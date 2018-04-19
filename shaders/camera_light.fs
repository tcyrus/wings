// $Id$
//
// Fragment shader for camera lighting
//
// Author: Dan Gudmundsson
//

#version 120

#ifdef GL_ARB_shader_texture_lod
#extension GL_ARB_shader_texture_lod : enable
#endif

uniform int UseDiffuseMap;
uniform int UseNormalMap;
uniform int UseMetallic;
uniform int UseEmission;
uniform int UseOcclusion;

uniform sampler2D DiffuseMap;
uniform sampler2D NormalMap;
uniform sampler2D RMMap;
uniform sampler2D EmissionMap;
uniform sampler2D OcculMap;

uniform sampler2D EnvBrdfMap;
uniform sampler2D EnvDiffMap;
uniform sampler2D EnvSpecMap;

uniform float metallic;
uniform float roughness;

varying vec3 ws_normal;
varying vec4 ws_tangent;
varying vec3 ws_position;
uniform vec3 ws_eyepoint;
varying vec4 color;
varying mat3 ws_TBN;

const vec3 l0_diff = vec3(0.7,0.7,0.7);
const vec3 l0_spec = vec3(0.2,0.2,0.2);
const vec3 l0_amb  = vec3(0.0,0.0,0.0);
const vec3 l0_pos  = vec3(0.0,0.5,0.0);

const vec3 lg_amb  = vec3(0.1,0.1,0.1);

const float M_PI = 3.141592653589793;
const float c_MinRoughness = 0.04;

// Encapsulate the various inputs used by the various functions in the shading equation
struct PBRInfo
{
    float NdotL;                  // cos angle between normal and light direction
    float NdotV;                  // cos angle between normal and view direction
    float NdotH;                  // cos angle between normal and half vector
    float LdotH;                  // cos angle between light direction and half vector
    float VdotH;                  // cos angle between view direction and half vector
    float perceptualRoughness;    // roughness value, as authored by the model creator (input to shader)
    float metalness;              // metallic value at the surface
    vec3 reflectance0;            // full reflectance color (normal incidence angle)
    vec3 reflectance90;           // reflectance color at grazing angle
    float alphaRoughness;         // roughness mapped to a more linear change in the roughness (proposed by [2])
    vec3 diffuseColor;            // color contribution from diffuse lighting
    vec3 specularColor;           // color contribution from specular lighting
};

vec4 SRGBtoLINEAR(vec4 srgbIn)
{
    // #ifdef MANUAL_SRGB
    // #ifdef SRGB_FAST_APPROXIMATION
    // vec3 linOut = pow(srgbIn.xyz,vec3(2.2));
    // #else //SRGB_FAST_APPROXIMATION
    vec3 bLess = step(vec3(0.04045),srgbIn.xyz);
    vec3 linOut = mix( srgbIn.xyz/vec3(12.92), pow((srgbIn.xyz+vec3(0.055))/vec3(1.055),vec3(2.4)), bLess );
    // #endif //SRGB_FAST_APPROXIMATION
    return vec4(linOut,srgbIn.w);;
    // #else //MANUAL_SRGB
    // return srgbIn;
    // #endif //MANUAL_SRGB
}

vec3 get_normal() {
    vec3 T = normalize(ws_tangent.xyz);
    float backface = (2.0 * float(gl_FrontFacing) - 1.0);
    if(UseNormalMap == 0 || dot(T,T) < 0.1)
	return backface*normalize(ws_normal); // No normal-map or Tangents
    // Calc Bumped normal
    vec3 BumpMapNormal = texture2D(NormalMap, gl_TexCoord[0].xy).xyz;
    BumpMapNormal = 2.0 * BumpMapNormal - vec3(1.0, 1.0, 1.0);
    vec3 N = normalize(ws_normal);
    T = normalize(T - dot(T, N) * N);
    vec3 B = cross(T, N) * ws_tangent.w;
    mat3 TBN = mat3(T, B, N);
    vec3 NewNormal = TBN * BumpMapNormal;
    // vec3 NewNormal = ws_TBN * BumpMapNormal;
    NewNormal = normalize(NewNormal);
    return backface*NewNormal;
}

vec4 get_diffuse() {
  if(UseDiffuseMap > 0) return SRGBtoLINEAR(texture2D(DiffuseMap, gl_TexCoord[0].xy));
  else return vec4(1.0, 1.0, 1.0, 1.0);
}

vec3 get_emission() {
  vec3 emi = vec3(gl_FrontMaterial.emission);
  if(UseEmission > 0) return SRGBtoLINEAR(texture2D(EmissionMap, gl_TexCoord[0].xy)).rgb * emi;
  return emi;
}

vec4 get_metalroughness() {
  vec4 mrSample = vec4(1.0,roughness,metallic,1.0);
  if(UseMetallic > 0) {
    mrSample *= texture2D(RMMap, gl_TexCoord[0].xy);
  }
  return clamp(mrSample, 0.0, 1.0);
}

// calculations

// The following equation models the Fresnel reflectance term of the spec equation (aka F())
// Implementation of fresnel from [4], Equation 15
vec3 specularReflection(PBRInfo pbrInputs)
{
    return pbrInputs.reflectance0 +
      (pbrInputs.reflectance90 - pbrInputs.reflectance0)
      * pow(clamp(1.0 - pbrInputs.VdotH, 0.0, 1.0), 5.0);
}

// This calculates the specular geometric attenuation (aka G()),
// where rougher material will reflect less light back to the viewer.
// This implementation is based on [1] Equation 4, and we adopt their modifications to
// alphaRoughness as input as originally proposed in [2].
float geometricOcclusion(PBRInfo pbrInputs)
{
    float NdotL = pbrInputs.NdotL;
    float NdotV = pbrInputs.NdotV;
    float r = pbrInputs.alphaRoughness;

    float attenuationL = 2.0 * NdotL / (NdotL + sqrt(r * r + (1.0 - r * r) * (NdotL * NdotL)));
    float attenuationV = 2.0 * NdotV / (NdotV + sqrt(r * r + (1.0 - r * r) * (NdotV * NdotV)));
    return attenuationL * attenuationV;
}

// The following equation(s) model the distribution of microfacet normals across the area being drawn (aka D())
// Implementation from "Average Irregularity Representation of a Roughened Surface for Ray Reflection" by T. S. Trowbridge, and K. P. Reitz
// Follows the distribution function recommended in the SIGGRAPH 2013 course notes from EPIC Games [1], Equation 3.
float microfacetDistribution(PBRInfo pbrInputs)
{
    float roughnessSq = pbrInputs.alphaRoughness * pbrInputs.alphaRoughness;
    float f = (pbrInputs.NdotH * roughnessSq - pbrInputs.NdotH) * pbrInputs.NdotH + 1.0;
    return roughnessSq / (M_PI * f * f);
}

// Basic Lambertian diffuse
// Implementation from Lambert's Photometria https://archive.org/details/lambertsphotome00lambgoog
// See also [1], Equation 1
vec3 diffuse(PBRInfo pbrInputs)
{
  return pbrInputs.diffuseColor / M_PI;
}

vec2 vec2uv(vec3 vec)
{
  float u = atan(vec.x, -vec.z);
  float v = -asin(vec.y);
  if(u < 0.0f) return vec2(0.5*u/M_PI, 0.5+v/M_PI);
  else return vec2(1.0f+0.5*u/M_PI, 0.5+v/M_PI);
}

vec3 background_ligthting(PBRInfo pbrInputs, vec3 N, vec3 reflection)
{
    vec3 refl = normalize(reflection);
    vec2 index = vec2uv(refl);

    //    if((bool)0)
    {
      // vec2 brdf = vec2(pbrInputs.NdotV, 1.0-pbrInputs.perceptualRoughness);
      // float difflight = (0.8*N.y+1.0)/2.0 * 0.15;
      // vec3 diffuse  = difflight * pbrInputs.diffuseColor;
      // float specligth = (0.8*reflection.y+1.0)/2.0 * 0.1;
      // vec3 specular = specligth * (pbrInputs.specularColor * brdf.x + brdf.y);
    // } else {
      vec2 brdf = texture2D(EnvBrdfMap, vec2(pbrInputs.NdotV, pbrInputs.perceptualRoughness)).rg;
      vec3 difflight = SRGBtoLINEAR(texture2D(EnvDiffMap, index)).rgb;
      vec3 diffuse  = difflight * pbrInputs.diffuseColor;
#ifdef GL_ARB_shader_texture_lod
      float mipCount = 8.0; // resolution of 512x256
      float lod = (pbrInputs.perceptualRoughness * mipCount);
      vec3 specligth = SRGBtoLINEAR(texture2DLod(EnvSpecMap, index, lod)).rgb;
#else
      vec3 specligth = SRGBtoLINEAR(texture2D(EnvSpecMap, index)).rgb;
#endif
      vec3 specular = specligth * (pbrInputs.specularColor * brdf.x + brdf.y);
      return diffuse + specular;
    }
}

void main(void)
{
    vec4 difftex = get_diffuse();
    vec4 baseColor = difftex*color;
    vec3 f0 = vec3(0.04);
    vec4 mr = get_metalroughness();
    float metallic = mr.b;
    float roughness = mr.g;
    // vec3 specular = gl_FrontMaterial.specular.rgb; // dgud
    vec3 specular = mix(f0, baseColor.rgb, metallic);

    vec3 n = get_normal();
    vec3 v = normalize(ws_eyepoint-ws_position);  // point to camera
    vec3 l = normalize(ws_eyepoint+l0_pos);  // point to ligth
    vec3 h = normalize(l+v);           // halfv between l and v

    float NdotL = clamp(dot(n, l), 0.001, 1.0);
    float NdotV = abs(dot(n, v)) + 0.001;
    float NdotH = clamp(dot(n, h), 0.0, 1.0);
    float LdotH = clamp(dot(l, h), 0.0, 1.0);
    float VdotH = clamp(dot(v, h), 0.0, 1.0);

    // Much inspired from glTF descriptions in: https://github.com/KhronosGroup/glTF-WebGL-PBR
    // http://blog.selfshadow.com/publications/s2012-shading-course/burley/s2012_pbs_disney_brdf_notes_v3.pdf
    float alphaRoughness = roughness*roughness;

    float reflectance = max(max(specular.r, specular.g), specular.b);
    float reflectance90 = clamp(reflectance * 25.0, 0.0, 1.0);
    vec3 specularEnvironmentR0 = specular.rgb;
    vec3 specularEnvironmentR90 = vec3(1.0, 1.0, 1.0) * reflectance90;

    PBRInfo pbrInputs = PBRInfo(
        NdotL,
        NdotV,
        NdotH,
        LdotH,
        VdotH,
        roughness,
        metallic,
        specularEnvironmentR0,
        specularEnvironmentR90,
        alphaRoughness,
	baseColor.rgb,
        specular
    );

    // Calculate the shading terms for the microfacet specular shading model
    vec3  F = specularReflection(pbrInputs);
    float G = geometricOcclusion(pbrInputs);
    float D = microfacetDistribution(pbrInputs);

    // Calculation of analytical lighting contribution
    vec3 diffuseContrib = (1.0 - max(max(F.r, F.g), F.b)) * diffuse(pbrInputs);
    vec3 specContrib = F * G * D / (4.0 * NdotL * NdotV);
    // Obtain final intensity as reflectance (BRDF) scaled by the energy of the light (cosine law)
    vec3 frag_color = NdotL * l0_diff * (diffuseContrib + specContrib);
    frag_color += 0.6*background_ligthting(pbrInputs, n, normalize(reflect(v, n)));

    if(UseOcclusion > 0) {
      float ao = texture2D(OcculMap, gl_TexCoord[0].xy).r;
      frag_color = mix(frag_color, frag_color * ao, 0.7);
    }

    frag_color += get_emission()*difftex.rgb;
    // frag_color = vec3(metallic);
    // frag_color = vec3(roughness);
    // frag_color = baseColor.rgb;
    // frag_color = diffuseContrib;
    // frag_color = specular;
    // frag_color = specContrib;
    // frag_color = F;
    // frag_color = 1.0 - F;
    // frag_color = 1.0 - vec3(max(max(F.r, F.g), F.b)),
    gl_FragColor = vec4(pow(frag_color, vec3(1.0/2.2)), baseColor.a); // Should be 2.2
    // frag_color = vec3(normalize(reflect(v, n)*+1.0)*0.5);
    // frag_color = n;
    // gl_FragColor = vec4((n+1.0)*0.5, baseColor.a);
}


