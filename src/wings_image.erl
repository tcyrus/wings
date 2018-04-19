%%
%%  wings_image.erl --
%%
%%     This module manages images.
%%
%%  Copyright (c) 2003-2011 Bjorn Gustavsson
%%
%%  See the file "license.terms" for information on usage and redistribution
%%  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
%%
%%     $Id$
%%

-module(wings_image).
-export([from_file/1,new/2,new_temp/2,new_hidden/2, create/1,
	 rename/2,txid/1,info/1,images/0,
	 screenshot/2,screenshot/1,viewport_screenshot/1,
	 bumpid/1,
	 is_normalmap/1,
	 next_id/0,delete_older/1,delete_from/1,delete/1,
	 update/2,update_filename/2,find_image/2,
	 window/1]).
-export([image_formats/0,image_read/1,image_write/1,
	 e3d_to_wxImage/1, wxImage_to_e3d/1]).
-export([maybe_exceds_opengl_caps/1]).

-behavior(gen_server).
-export([start_link/0, init/1, handle_call/3, handle_cast/2]).

-define(NEED_OPENGL, 1).
-include("wings.hrl").
-include("e3d_image.hrl").
-import(lists, [reverse/1]).

start_link() ->
    Env = wings_io:get_process_option(),
    gen_server:start_link({local, ?MODULE}, ?MODULE, [Env], [{spawn_opt, [{fullsweep_after,0}]}]).

%%%
%%% Interface against plug-ins.
%%%
image_formats() ->
    wings_plugin:call_ui({image,formats,[]}).

image_read(Ps) ->
    CurrDir  = wings_pref:get_value(current_directory),
    OptDir   = proplists:get_value(opt_dir,  Ps, undefined),
    AbsFile  = proplists:get_value(filename, Ps),
    FileName = filename:basename(AbsFile),
    Dirs     = filename:split(filename:dirname(AbsFile)),
    case wings_plugin:call_ui({image,read,Ps}) of
	{error, enoent} ->
	    try_relative_paths(Dirs,FileName,OptDir,CurrDir,Ps);
	Ok = #e3d_image{} ->
	    Ok#e3d_image{filename=AbsFile};
	Other ->
	    Other
    end.

image_write(Ps) ->
    case catch wings_plugin:call_ui({image,write,Ps}) of
	{'EXIT',Reason} ->
	    {error,{none,?MODULE,{crash,Reason}}};
	Result -> Result
    end.

e3d_to_wxImage(I) ->
    e3d_to_wxImage_1(I).

%%%
%%% Client API.
%%%

from_file(Filename) ->
    Props = [{filename,Filename},{alignment,1}],
    case image_read(Props) of
	#e3d_image{}=Image ->
	    Name = filename:basename(Filename),
	    req({new,Image#e3d_image{name=Name},false, false});
	{error,_}=Error -> Error
    end.

new(Name, E3DImage) ->
    req({new,E3DImage#e3d_image{name=Name},false, false}).

new_temp(Name, E3DImage) ->
    req({new,E3DImage#e3d_image{name=Name},true, false}).

new_hidden(Name, E3DImage) ->
    req({new,E3DImage#e3d_image{name=Name},true, true}).

create(St) ->
    create_image(),
    St.

rename(Id, NewName) ->
    req({rename,Id,NewName}).

txid(Id) ->
    req({txid,Id}, false).

is_normalmap(Id) ->
    req({is_normalmap,Id},false).

bumpid(Id) ->
    req({bumpid,Id}, false).

info(Id) ->
    req({info,Id}, false).

images() ->
    req(images, false).

next_id() ->
    req(next_id, false).

delete_older(Id) ->
    req({delete_older,Id}).

delete_from(Id) ->
    req({delete_from,Id}).

delete(Id) ->
    req({delete,Id}).

update(Id, Image) ->
    req({update,Id,Image}).

update_filename(Id, Filename) ->
    gen_server:cast(?MODULE, {update_filename,Id,Filename}).

find_image(Dir, Filename) ->
    req({find_image,Dir,Filename}, false).

req(Req) ->
    req(Req, true).

req(Req, Notify) ->
    Reply = gen_server:call(?MODULE, Req),
    Running = get(wings_not_running) == undefined,
    case Notify andalso Running of
        false -> ok;
        true -> wings_wm:notify(image_change)
    end,
    Reply.

screenshot(Ask, _) when is_atom(Ask) ->
    ViewPortOnly = wings_pref:get_value(screenshot_viewport_only, false),
    SaveView = wings_pref:get_value(screenshot_save_current_view, false),
    Qs = [{?__(2,"Capture viewport only"),ViewPortOnly},
          {?__(3,"Add current view to Saved Views"),SaveView},
          {hframe,[{label,?__(4,"Name")},
		   {text,?__(1,"Screenshot"),[]}]}],
    wings_dialog:dialog(Ask, ?__(1,"Screenshot"), [{vframe,Qs}],
			fun(Res) -> {tools,{screenshot,Res}} end);
screenshot([ViewPortOnly,SaveView,Name], St) ->
    wings_pref:set_value(screenshot_viewport_only, ViewPortOnly),
    wings_pref:set_value(screenshot_save_current_view, SaveView),
    wings_wm:send_after_redraw(geom,{action,{tools,{screenshot,[ViewPortOnly,Name]}}}),
    case SaveView of
      true -> wings_view:command({views,{save,[Name]}},St);
      false -> St
    end;
screenshot([ViewPortOnly,Name], St) ->
    case ViewPortOnly of
      true -> viewport_screenshot(Name);
      false -> screenshot(Name)
    end,
    St.

screenshot(Name) ->
    {W,H} = wings_wm:top_size(),
    gl:pixelStorei(?GL_PACK_ALIGNMENT, 1),
    gl:readBuffer(?GL_FRONT),
    Mem = wings_io:get_buffer(W*H*3, ?GL_UNSIGNED_BYTE),
    gl:readPixels(0, 0, W, H, ?GL_RGB, ?GL_UNSIGNED_BYTE, Mem),
    ImageBin = wings_io:get_bin(Mem),
    Image = #e3d_image{image=ImageBin,width=W,height=H},
    Id = new_temp(Name, Image),
    window(Id).

viewport_screenshot(Name) ->
%% screenshot of just the viewport scene
    {X,Y,W,H} = wings_wm:viewport(),
    gl:pixelStorei(?GL_PACK_ALIGNMENT, 1),
    gl:readBuffer(?GL_FRONT),
    Mem = wings_io:get_buffer(W*H*3, ?GL_UNSIGNED_BYTE),
    gl:readPixels(X, Y, W, H, ?GL_RGB, ?GL_UNSIGNED_BYTE, Mem),
    ImageBin = wings_io:get_bin(Mem),
    Image = #e3d_image{image=ImageBin,width=W,height=H},
    Id = new_temp(Name, Image),
    window(Id).

%%%
%%% Server implementation.
%%%
%%% Reason for using a server: Convenient encapsulation of images.
%%% We'll get an entire process dictionary to use for texture ids.
%%%

-record(ist,
	{next=0,				%Next image ID.
	 images					%All images (gb_trees).
	}).

init([Env]) ->
    process_flag(trap_exit, true),
    wings_io:set_process_option(Env),
    {ok, #ist{images=gb_trees:empty()}}.

handle_call({new,#e3d_image{name=Name0}=Im0,false,Hide}, _From, #ist{next=Id,images=Images0}=S) ->
    Name = make_unique(Name0, Images0),
    Im = maybe_convert(Im0#e3d_image{name=Name}),
    Images = gb_trees:insert(Id, hide(Im, Hide), Images0),
    case make_texture(Id, Im) of
	{error,_GlErr}=Err ->
	    {reply, Err, S};
	_ ->
	    {reply, Id, S#ist{next=Id+1,images=Images}}
    end;
handle_call({new,#e3d_image{name=Name}=Im,true,Hide}, From, #ist{images=Images}=S0) ->
    Exist = fun({_, #e3d_image{name=N}}) when N =:= Name -> true;
	       ({_, {hidden, #e3d_image{name=N}}}) when N =:= Name -> true;
	       (_) -> false
	    end,
    case lists:filter(Exist, gb_trees:to_list(Images)) of
	[] ->
	    handle_call({new,Im,false,Hide}, From, S0);
	[{Id,_}] ->
	    {reply, ok, S} = handle_call({delete,Id}, From, S0),
	    handle_call({new,Im,false,Hide}, From, S)
    end;
handle_call({rename,Id,Name0}, _From, #ist{images=Images0}=S) ->
    Name = make_unique(Name0, gb_trees:delete(Id, Images0)),
    Im0 = gb_trees:get(Id, Images0),
    Im = (image_rec(Im0))#e3d_image{name=Name},
    Images = gb_trees:update(Id, hide(Im, is_hidden(Im0)), Images0),
    {reply, Id,S#ist{images=Images}};
handle_call({txid,Id}, _From, S) ->
    case get(Id) of
        undefined -> {reply, none, S};
        TxId -> {reply, TxId, S}
    end;
handle_call({bumpid,Id}, _From, S) ->
    case get({Id,bump}) of
        undefined -> {reply, create_bump(Id,undefined,S), S};
        TxId -> {reply, TxId, S}
    end;
handle_call({is_normalmap,Id}, _From, S) ->
    case {get({Id,bump}),get(Id)} of
        {undefined,undefined} -> {reply, none, S};
        {undefined, ImId} -> {reply, create_bump(Id,ImId,S), S};
        {TxId,_} -> {reply, TxId, S}
    end;
handle_call({info,Id}, _From, #ist{images=Images}=S) ->
    case gb_trees:lookup(Id, Images) of
	{value,E3D} -> {reply, image_rec(E3D),S};
	none -> {reply, none,S}
    end;
handle_call(images, _From, #ist{images=Images}=S) ->
    {reply, [NotHidden || NotHidden = {_, #e3d_image{}} <- gb_trees:to_list(Images)],S};
handle_call(next_id, _From, #ist{next=Id}=S) ->
    {reply,Id,S};
handle_call({delete,Id}, _From, #ist{images=Images0}=S) ->
    delete_bump(Id),
    wings_gl:deleteTextures([erase(Id)]),
    Images = gb_trees:delete(Id, Images0),
    {reply, ok, S#ist{images=Images}};
handle_call({delete_older,Id}, _From, #ist{images=Images0}=S) ->
    Images1 = delete_older_1(gb_trees:to_list(Images0), Id),
    Images = gb_trees:from_orddict(Images1),
    {reply, ok, S#ist{images=Images}};
handle_call({delete_from,Id}, _From, #ist{images=Images0}=S) ->
    Images1 = delete_from_1(gb_trees:to_list(Images0), Id, []),
    Images = gb_trees:from_orddict(Images1),
    {reply, ok, S#ist{images=Images}};
handle_call({update,Id,Image}, _From, S) ->
    do_update(Id, Image, S);
handle_call({find_image, Dir, File}, _From, #ist{images=Ims}=S) ->
    AbsName = filename:join(Dir, File),
    Find = fun(Fn) -> Fn == AbsName end,
    Found = [Id || {Id, #e3d_image{filename=FN}} <-
                       gb_trees:to_list(Ims), Find(FN)],
    case Found of
        [] -> {reply, false, S};
        [Id|_] -> {reply, {true, Id}, S}
    end;

%% handle_call({load_envmap, Img}, _From, #ist{} = S) ->
%%     case cl_setup() of
%%         {error, _R} = Err -> {reply, Err, S};
%%         CL  -> {reply, make_envmap(CL, Img), S}
%%     end;
handle_call(Req, _From, S) ->
    io:format("~w: Bad request: ~w~n", [?MODULE, Req]),
    {reply, error, S}.

handle_cast({update_filename,Id,NewName}, #ist{images=Images0}=S) ->
    Im0 = gb_trees:get(Id, Images0),
    Im = (image_rec(Im0))#e3d_image{filename=NewName},
    Images = gb_trees:update(Id, hide(Im, is_hidden(Im0)), Images0),
    {noreply, S#ist{images=Images}}.


%%%%%%%%%%%%%%%%% Internal Functions %%%%%%%%%%%%%%%


create_bump(Id, BumpId, #ist{images=Images0}) ->
    delete_bump(Id),  %% update case..
    case gb_trees:lookup(Id, Images0) of
	{value, E3D0} ->
	    E3D = image_rec(E3D0),
	    gl:pushAttrib(?GL_TEXTURE_BIT),
	    case get(Id) of
		BumpId ->
		    gl:bindTexture(?GL_TEXTURE_2D, BumpId),
		    Image = case e3d_image:convert(maybe_scale(E3D),r8g8b8,1,lower_left) of
				E3D -> E3D;
				New = #e3d_image{width=W,height=H,image=Bits} ->
				    gl:texImage2D(?GL_TEXTURE_2D,0,?GL_RGB,W,H,0,?GL_RGB,
						  ?GL_UNSIGNED_BYTE,Bits),
				    New
			    end,
		    MipMapFun = fun() -> e3d_image:buildNormalMipmaps(Image) end,
		    MipMaps = worker_process(MipMapFun),
		    TxId = BumpId;
		_ ->
		    %% Scale ?? 4 is used in the only example I've seen.
		    Img = e3d_image:convert(maybe_scale(E3D), r8g8b8, 1, lower_left),
		    {#e3d_image{width=W,height=H,image=Bits},MipMaps}
			= e3d_image:height2normal(Img, #{scale=>4.0}, true),
		    [TxId] = gl:genTextures(1),
		    gl:bindTexture(?GL_TEXTURE_2D, TxId),
		    gl:texImage2D(?GL_TEXTURE_2D,0,?GL_RGB,W,H,0,?GL_RGB,
				  ?GL_UNSIGNED_BYTE,Bits)
	    end,
	    put({Id,bump}, TxId),
	    update_mipmaps(TxId,MipMaps),
	    gl:popAttrib(),
	    TxId;
	_ ->
	    none
    end.

update_mipmaps(TxId, MipMaps) ->
    gl:enable(?GL_TEXTURE_2D),
    gl:bindTexture(?GL_TEXTURE_2D, TxId),
    gl:texParameteri(?GL_TEXTURE_2D, ?GL_TEXTURE_MAG_FILTER, ?GL_LINEAR),
    gl:texParameteri(?GL_TEXTURE_2D, ?GL_TEXTURE_MIN_FILTER, ?GL_LINEAR_MIPMAP_LINEAR),
    gl:texParameteri(?GL_TEXTURE_2D, ?GL_TEXTURE_WRAP_S, ?GL_REPEAT),
    gl:texParameteri(?GL_TEXTURE_2D, ?GL_TEXTURE_WRAP_T, ?GL_REPEAT),
    Load = fun({Level,MW,MH, Bin}) ->
		   gl:texImage2D(?GL_TEXTURE_2D,Level,?GL_RGB,MW,MH,0,
				 ?GL_RGB, ?GL_UNSIGNED_BYTE, Bin)
	   end,
    [Load(MM) || MM <- MipMaps].


maybe_convert(#e3d_image{type=Type0,order=Order}=Im) ->
    case {img_type(Type0),Order} of
	{Type0,lower_left} -> Im;
	{Type,_} -> e3d_image:convert(Im, Type, 1, lower_left)
    end.

img_type(b8g8r8) -> r8g8b8;
img_type(b8g8r8a8) -> r8g8b8a8;
img_type(Type) -> Type.

make_texture(Id, Image) ->
    case init_texture(Image) of
	{error,_}=Error ->
	    Error;
	TxId ->
	    put(Id, TxId),
	    TxId
    end.

init_texture(Image) ->
    case get(wings_not_running) of
        true -> 0;
        _ ->
            [TxId] = gl:genTextures(1),
            case init_texture(image_rec(Image), TxId) of
                {error,_}=Error ->
                    wings_gl:deleteTextures([TxId]),
                    Error;
                Other ->
                    Other
            end
    end.

init_texture(Image0, TxId) ->
    case maybe_scale(Image0) of
        {error,_}=Error ->
            Error;
        Image ->
            #e3d_image{width=W,height=H,image=Bits} = Image,
            gl:pushAttrib(?GL_TEXTURE_BIT),
            gl:enable(?GL_TEXTURE_2D),
            gl:bindTexture(?GL_TEXTURE_2D, TxId),
            Ft=case wings_pref:get_value(filter_texture, false) of
				true -> ?GL_LINEAR;
				false -> ?GL_NEAREST
            end,
            case wings_gl:is_ext({1,4},'GL_SGIS_generate_mipmap') of
		true ->
		    gl:texParameteri(?GL_TEXTURE_2D, ?GL_GENERATE_MIPMAP, ?GL_TRUE),
		    gl:texParameteri(?GL_TEXTURE_2D, ?GL_TEXTURE_MIN_FILTER,
				     ?GL_LINEAR_MIPMAP_LINEAR);
		false ->
		    gl:texParameteri(?GL_TEXTURE_2D, ?GL_TEXTURE_MIN_FILTER, Ft)
            end,
            gl:texParameteri(?GL_TEXTURE_2D, ?GL_TEXTURE_MAG_FILTER, Ft),
            gl:texParameteri(?GL_TEXTURE_2D, ?GL_TEXTURE_WRAP_S, ?GL_REPEAT),
            gl:texParameteri(?GL_TEXTURE_2D, ?GL_TEXTURE_WRAP_T, ?GL_REPEAT),
            Format = texture_format(Image),
            gl:texImage2D(?GL_TEXTURE_2D, 0, internal_format(Format),
			  W, H, 0, Format, ?GL_UNSIGNED_BYTE, Bits),
            gl:popAttrib(),
            TxId
    end.

maybe_scale(#e3d_image{width=W0,height=H0}=Image) ->
%%  case wings_gl:is_ext({2,0}, 'GL_ARB_texture_non_power_of_two') of
%%  Aarg ATI doesn't support ARB_NPOT textures, though it report GL_VER >= 2.0
    case maybe_exceds_opengl_caps(Image) of
        {error,_}=Error ->
            Error;
        Image1 ->
            #e3d_image{width=W1,height=H1}=Image1,
            case {W1,H1} of
                {W0,H0} ->
                    GL_ARB = wings_gl:is_ext('GL_ARB_texture_non_power_of_two');
                {_,_} -> GL_ARB = false
            end,
            case GL_ARB of
                true ->
                    Image;
                false ->
                    case {nearest_power_two(W1),nearest_power_two(H1)} of
                        {W1,H1} ->
                            Image1;
                        {W,H} ->
                            resize_image(Image1, W, H)
                    end
            end
    end.

maybe_exceds_opengl_caps(#e3d_image{width=W0,height=H0}=Image) ->
    MaxSize = hd(gl:getIntegerv(?GL_MAX_TEXTURE_SIZE)),
    case need_resize_image(W0, H0, MaxSize) of
        true ->
            ScaleFactor = case W0 > H0 of
                true ->
                    MaxSize/W0;
                false ->
                    MaxSize/H0
            end,
            W = trunc(W0*ScaleFactor),
            H = trunc(H0*ScaleFactor),
            resize_image(Image, W, H);
        false ->
            Image
    end.

resize_image(#e3d_image{width=W0,height=H0,bytes_pp=BytesPerPixel,
  image=Bits0}=Image, W, H) ->
    Out = wings_io:get_buffer(BytesPerPixel*W*H, ?GL_UNSIGNED_BYTE),
    Format = texture_format(Image),
    GlErr =glu:scaleImage(Format, W0, H0, ?GL_UNSIGNED_BYTE,
        Bits0, W, H, ?GL_UNSIGNED_BYTE, Out),
    case GlErr of
        0 ->
            Bits = wings_io:get_bin(Out),
            Image#e3d_image{width=W,height=H,bytes_pp=BytesPerPixel,image=Bits};
        _ ->
            {error,GlErr}
    end.

need_resize_image(W, H, Max) when W > Max; H > Max ->
    true;
need_resize_image(_, _, _) ->
    false.

nearest_power_two(N) when (N band -N) =:= N -> N;
nearest_power_two(N) -> nearest_power_two(N, 1).

nearest_power_two(N, B) when B > N -> B bsr 1;
nearest_power_two(N, B) -> nearest_power_two(N, B bsl 1).

texture_format(#e3d_image{type=r8g8b8}) -> ?GL_RGB;
texture_format(#e3d_image{type=r8g8b8a8}) -> ?GL_RGBA;
texture_format(#e3d_image{type=b8g8r8}) -> ?GL_BGR;
texture_format(#e3d_image{type=b8g8r8a8}) -> ?GL_BGRA;
texture_format(#e3d_image{type=g8}) -> ?GL_LUMINANCE;
texture_format(#e3d_image{type=a8}) -> ?GL_ALPHA.

%% Long ago we tried to use compression, but we no longer do since compression
%% lowers the quality too much, especially for bump/normal maps.
internal_format(?GL_BGR) -> ?GL_RGB;
internal_format(?GL_BGRA) -> ?GL_RGBA;
internal_format(Else) -> Else.

delete_older_1([{Id,_}|T], Limit) when Id < Limit ->
    delete_bump(Id),
    wings_gl:deleteTextures([erase(Id)]),
    delete_older_1(T, Limit);
delete_older_1(Images, _) -> Images.

delete_from_1([{Id,_}=Im|T], Limit, Acc) when Id < Limit ->
    delete_from_1(T, Limit, [Im|Acc]);
delete_from_1([{Id,_}|T], Limit, Acc) ->
    delete_bump(Id),
    wings_gl:deleteTextures([erase(Id)]),
    delete_from_1(T, Limit, Acc);
delete_from_1([], _, Acc) -> reverse(Acc).

delete_bump(Id) ->
    TxId = get(Id), 
    case erase({Id,bump}) of
	undefined -> ok;
	TxId -> ok;
	Bid ->  wings_gl:deleteTextures([Bid])
    end.

do_update(Id, In = #e3d_image{width=W,height=H,type=Type,name=NewName}, 
	  #ist{images=Images0}=S) ->
    Image = gb_trees:get(Id, Images0),
    Im0 = #e3d_image{filename=File,name=OldName} = image_rec(Image),
    Name = if is_list(NewName), length(NewName) > 2 -> NewName;
	      true -> OldName
	   end,
    Im   = maybe_convert(In#e3d_image{filename=File, name=Name}),
    TxId = get(Id),
    Images = gb_trees:update(Id, hide(Im, is_hidden(Image)), Images0),
    Size = {Im0#e3d_image.width, Im0#e3d_image.height, Im0#e3d_image.type},
    case Size of
	{W,H,Type} ->
	    gl:bindTexture(?GL_TEXTURE_2D, TxId),
	    gl:texSubImage2D(?GL_TEXTURE_2D, 0, 0, 0,
			     W, H, texture_format(Im), 
			     ?GL_UNSIGNED_BYTE, Im#e3d_image.image);
	_ ->
	    init_texture(Im, TxId)
    end,
    case get({Id,bump}) of
	undefined ->
	    S#ist{images=Images};
	Bid ->
	    create_bump(Id, Bid, S#ist{images=Images}),
	    S#ist{images=Images}
    end.

make_unique(Name, Images0) ->
    Images = [(image_rec(Image))#e3d_image.name || Image <- gb_trees:values(Images0)],
    wings_util:unique_name(Name, Images).

hide(Im, true) -> {hidden, Im};
hide(Im, false) -> Im.

image_rec({hidden, Im}) -> Im;
image_rec(#e3d_image{} = Im) -> Im.

is_hidden({hidden, _}) -> true;
is_hidden(#e3d_image{}) -> false.

%%%
%%% Window for image.
%%%

window(Id) ->
    Name = {image,Id},
    case wings_wm:is_window(Name) of
	true ->
	    wings_wm:raise(Name);
	false ->
	    wings_image_viewer:new(Name, info(Id)),
	    keep
    end.

%%%
%%% Creating images with pre-defined patterns.
%%%

create_image() ->
    Def = wings_pref:get_value(current_directory),
    Ps = [{extensions,image_formats()},{multiple,false}],
    Flags = [{props, Ps}],
    Qs = [{label_column,
           [{?__(0,"Import"), {button, {text, Def, Flags}}},
            separator,
            {?__(1,"Width"), {text, 256,[{range,{8,1024}}]}},
            {?__(2,"Height"),{text, 256,[{range,{8,1024}}]}}
           ]},
	  {vradio,
	   [{?__(4,"Grid"),grid},
	    {?__(5,"Checkerboard"),checkerboard},
	    {?__(6,"Vertical Bars"),vbars},
	    {?__(7,"Horizontal Bars"),hbars},
	    {?__(8,"White"),white},
	    {?__(9,"Black"),black}],
	   grid, [{title, ?__(3,"Pattern")}]}],
    wings_dialog:dialog(?__(10,"Create Image"), Qs,
			fun([File, W,H,Pattern]) ->
                                case filelib:is_regular(File) of
                                    true  ->
                                        from_file(File);
                                    false ->
                                        create_image_1(Pattern, W, H)
                                end,
				ignore
			end).

create_image_1(Pattern, W, H) ->
    Pixels = pattern(Pattern, W, H),
    case {byte_size(Pixels),3*W*H} of
	{S,S} ->				%Assertion.
	    Im = #e3d_image{width=W,height=H,image=Pixels,order=upper_left},
	    new(atom_to_list(Pattern), Im)
    end.

pattern(grid, W, H) ->
    grid(W, H);
pattern(checkerboard, W, H) ->
    checkerboard(W, H);
pattern(vbars, W, H) ->
    vertical_bars(W, H);
pattern(hbars, W, H) ->
    horizontal_bars(W, H);
pattern(white, W, H) ->
    all_white(W, H);
pattern(black, W, H) ->
    all_black(W, H).

%% Generate a grid image.
grid(Width, Height) ->
    White = [255,255,255],
    Black = [0,0,0],
    WhiteRow = pattern_repeat(Width, White),
    BlackLine = pattern_repeat(14, Black),
    Mixed0 = pattern_repeat(14*((Width+15) div 16), [White,BlackLine|White]),
    Mixed = truncate(Mixed0, 3*Width),
    MixedRows = pattern_repeat(14, Mixed),
    R = [WhiteRow,MixedRows|WhiteRow],
    All = pattern_repeat((Height+15) div 16, R),
    truncate(All, 3*Width*Height).

%% Generate a checkerboard image of 4x4 squares 
%% with given side length in pixels.
checkerboard(Width, Height) ->
    FourWhite = pattern_repeat(3*4, 255),
    FourBlack = pattern_repeat(3*4, 0),
    RoundedW = (Width+7) div 8,
    RowSize = 3*Width,
    R1 = truncate(pattern_repeat(RoundedW, [FourBlack|FourWhite]), RowSize),
    R2 = truncate(pattern_repeat(RoundedW, [FourWhite|FourBlack]), RowSize),
    R8 = [pattern_repeat(4, [R1])|pattern_repeat(4, [R2])],
    truncate(pattern_repeat(RoundedW, R8), 3*Width*Height).

%% Generate a vertical bars image of 4 pixels width of given size.
vertical_bars(Width, Height) ->
    W4 = pattern_repeat(3*4, 255),
    B4 = pattern_repeat(3*4, 0),
    Row = truncate(pattern_repeat((Width+7) div 8, [B4|W4]), 3*Width),
    list_to_binary(pattern_repeat(Height, Row)).

%% Generate a horizontal bars image of 4 pixels height.
horizontal_bars(Width, Height) ->
    WhiteRow = pattern_repeat(3*Width, 255),
    BlackRow = pattern_repeat(3*Width, 0),
    WR4 = pattern_repeat(4, WhiteRow),
    BR4 = pattern_repeat(4, BlackRow),
    truncate(pattern_repeat((Height+7) div 8, [BR4|WR4]), 3*Width*Height).

%% Generate an all white image with given size.
all_white(Width, Height) ->
    solid(Width, Height, 255).

%% Generate an all white image with given size.
all_black(Width, Height) ->
    solid(Width, Height, 0).

solid(Width, Height, Channel) ->
    list_to_binary(pattern_repeat(3*Width*Height, Channel)).

pattern_repeat(0, _) -> [];
pattern_repeat(1, D) -> [D];
pattern_repeat(N, D) ->
    B = pattern_repeat(N div 2, D),
    case N rem 2 of
	0 -> [B|B];
	1 -> [D,B|B]
    end.

truncate(B0, Sz) ->
    <<B:Sz/binary,_/binary>> = list_to_binary(B0),
    B.

%% Run a computation in a worker process with a generous heap size
%% and the default generational garbage collector. Before R12B,
%% heap fragments would allow the heap to be over-committed, but in
%% R12B there will be garbage collection as soon as the heap space
%% is exhausted.
worker_process(WorkFun) ->
    ResultRef = make_ref(),
    {Pid,Ref} = spawn_opt(fun() -> exit({ResultRef,WorkFun()}) end,
			  [monitor,{min_heap_size,20000}]),
    receive
	{'DOWN',Ref,process,Pid,{ResultRef,Res}} ->
	    Res;
	{'DOWN',Ref,process,Pid,Reason} ->
	    exit(Reason)
    end.

e3d_to_wxImage_1(I = #e3d_image{bytes_pp=4, width=W, height=H}) ->
    #e3d_image{image=RGB} = e3d_image:convert(I, r8g8b8, 1, upper_left),
    Wx = wxImage:new(W,H,RGB),
    #e3d_image{image=Alpha} = e3d_image:convert(I, a8, 1, upper_left),
    wxImage:setAlpha(Wx, Alpha),
    Wx;
e3d_to_wxImage_1(I = #e3d_image{bytes_pp=3, width=W, height=H}) ->
    #e3d_image{image=RGB} = e3d_image:convert(I, r8g8b8, 1, upper_left),
    wxImage:new(W,H,RGB);
e3d_to_wxImage_1(I = #e3d_image{bytes_pp=1, width=W, height=H}) ->
    #e3d_image{image=RGB} = e3d_image:convert(I, r8g8b8, 1, upper_left),
    wxImage:new(W,H,RGB).

wxImage_to_e3d(Wx) ->
    wxImage_to_e3d_1(Wx).

wxImage_to_e3d_1(Wx) ->
    wxImage = wx:getObjectType(Wx), %% Assert
    E3d0 = #e3d_image{image=wxImage:getData(Wx),
		      width=wxImage:getWidth(Wx),
		      height=wxImage:getHeight(Wx),
		      order = upper_left
		     },
    case wxImage:hasAlpha(Wx) of
	true ->
            e3d_image:add_alpha(E3d0, wxImage:getAlpha(Wx));
	false ->
            case wxImage:hasMask(Wx) of
                true ->
                    wxImage:initAlpha(Wx),
                    e3d_image:add_alpha(E3d0, wxImage:getAlpha(Wx));
                false ->
                    E3d0
            end
    end.

try_relative_paths(Dirs,FileName,undefined,CurrDir,Ps) ->
    try_relative_paths(CurrDir,Dirs,FileName,Ps);
try_relative_paths(Dirs,FileName,OptDir,CurrDir,Ps) ->
    case try_relative_paths(OptDir,Dirs,FileName,Ps) of
	{error, enoent} ->
	    try_relative_paths(Dirs,FileName,undefined,CurrDir,Ps);
	Other ->
	    Other
    end.

try_relative_paths(Start, Rel, File, Ps0) ->
    Abs = case Rel of
	      [] -> filename:join(Start,File);
	      _  -> filename:join(filename:join(Start,filename:join(Rel)),File)
	  end,
    Ps = lists:keyreplace(filename,1,Ps0,{filename,Abs}),
    case wings_plugin:call_ui({image,read,Ps}) of
	Err = {error, enoent} ->
	    case Rel of
		[] -> Err;
		[_|Rest] -> try_relative_paths(Start,Rest,File,Ps0)
	    end;
	Ok = #e3d_image{} ->
	    Ok#e3d_image{filename=Abs};
	Other ->
	    Other
    end.
