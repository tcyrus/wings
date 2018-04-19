// $Id$
//
// Vertex shader for hemispherical lighting
//
// Author: Randi Rost
//
// Copyright (C) 2005 3Dlabs, Inc.
//
// See 3Dlabs-License.txt for license information
//

attribute vec4 wings_tangent;

varying vec3 normal;      // eye space
varying vec3 ws_normal;   // world space
varying vec4 ws_tangent;  // world space
varying vec3 ws_position; // world space
varying vec3 ecPosition;  // eye space
varying vec4 color;
varying mat3 ws_TBN;      // world space

void main(void)
{
    ws_position = gl_Vertex;
    ecPosition  = vec3(gl_ModelViewMatrix * gl_Vertex);
    color	= gl_FrontMaterial.diffuse * gl_Color;
    normal	= normalize(gl_Normal);
    vec3 T      = normalize(wings_tangent.xyz);
    // T = normalize(T - dot(T, normal) * normal); // ??
    vec3 B = cross(T, normal) * wings_tangent.w;
    ws_TBN = mat3(T, B, normal);
    ws_normal = normal;
    ws_tangent = vec4(T, wings_tangent.w);
    normal = gl_NormalMatrix * normal;
    gl_TexCoord[0]	= gl_MultiTexCoord0;
#ifdef __GLSL_CG_DATA_TYPES // Fix clipping for Nvidia and ATI
    gl_ClipVertex   = gl_ModelViewMatrix * gl_Vertex;
#endif
    gl_Position 	= ftransform();
}
