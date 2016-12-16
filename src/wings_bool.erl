%%
%%  wings_bool.erl --
%%
%%     This module implements boolean commands for Wings objects.
%%
%%  Copyright (c) 2016 Dan Gudmundsson
%%
%%  See the file "license.terms" for information on usage and redistribution
%%  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
%%
%%     $Id$
%%

-module(wings_bool).
-export([add/1]).

-include("wings.hrl").
-compile(export_all).
-define(EPSILON, 1.0e-10).  %% used without SQRT() => 1.0e-5

add(#st{shapes=Sh0}=St0) ->
    MapBvh = wings_sel:fold(fun make_bvh/3, [], St0),
    try
        IntScts = find_intersects(MapBvh, MapBvh, []),
        Sel = [{Id,gb_sets:from_list(Es)} || #{id:=Id, es:=Es} <- IntScts],
        %% ?dbg("~p~n", [Sel2]),
        Upd = fun(#{we:=#we{id=Id}=We, delete:=Del}, Sh) ->
                      gb_trees:update(Id, We, gb_trees:delete(Del, Sh))
              end,
        Sh = lists:foldl(Upd, Sh0, IntScts),
        wings_sel:set(edge, Sel, St0#st{shapes=Sh})
    catch throw:multiple_face_ints ->
            Str = ?__(1, "Have not yet learned how to handle multiple intersections"
                      " in the same face,\ntry subdividing intersected faces first."),
            wings_u:error_msg(Str)
    end.

find_intersects([Head|Tail], [_|_] = Alts, Acc) ->
    case find_intersect(Head, Alts) of
	false ->
	    find_intersects(Tail,Alts, Acc);
	{H1,Merge} ->
	    NewTail = remove(H1, Head, Tail),
	    NewAlts = remove(H1, Head, Alts),
	    find_intersects(NewTail, NewAlts, [Merge|remove(Head,Acc)])
    end;
find_intersects(_, _, Acc) -> lists:reverse(Acc).

find_intersect(#{id:=Id}=Head, [#{id:=Id}|Rest]) ->
    find_intersect(Head, Rest);
find_intersect(#{bvh:=B1}=Head, [#{bvh:=B2}=H1|Rest]) ->
    case e3d_bvh:intersect(B1, B2) of
	[] ->  find_intersect(Head,Rest);
	EdgeInfo -> {H1, merge(EdgeInfo, Head, H1)}
    end;
find_intersect(_Head, []) ->
    false.

merge(EdgeInfo, #{id:=Id1,we:=We1}=I1, #{id:=Id2,we:=We2}=I2) ->
    %?dbg("~p~n",[EdgeInfo]),
    ReEI0 = [remap(Edge, I1, I2) || Edge <- EdgeInfo],
    %?dbg("~p~n",[ReEI0]),
    {Vmap, ReEI} = make_vmap(ReEI0, We1, We2),
    Loops0 = build_vtx_loops(ReEI, []),
    %?dbg("~p~n",[Loops0]),
    L10 = [split_loop(Loop, Vmap, {We1,We2}) || Loop <- Loops0],
    L20 = [split_loop(Loop, Vmap, {We2,We1}) || Loop <- Loops0],
    Loops = [filter_tri_edges(Loop,We1,We2) || Loop <- lists:zip(L10,L20)],

    {L1,L2} = lists:unzip(Loops),

    {Es1, We1N0} = make_verts(L1, Vmap, We1),
    {Es2, We2N0} = make_verts(L2, Vmap, We2),
    %% ?dbg("loops~n",[]), [io:format("  ~w~n",[Fs]) || {_, Fs} <- Es2],
    DRes1 = dissolve_faces_in_edgeloops(Es1, We1N0),
    DRes2 = dissolve_faces_in_edgeloops(Es2, We2N0),

    {We,Es} = weld(DRes1, DRes2),
    [Del] = lists:delete(We#we.id, [Id1,Id2]),
    ok = wings_we_util:validate(We),
    #{id=>We#we.id, es=>Es, we=>We, delete=>Del}.
    %#{id=>Id2, es=>Es2, we=>We2N0, delete=>Id1}. %% debug

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dissolve_faces_in_edgeloops(Es, #we{fs=_Ftab} = We0) ->
    Fs = select_region(Es, We0),
    We = wings_dissolve:faces(Fs, We0),
    Faces = wings_we:new_items_as_ordset(face, We0, We),
    {order_loops(Faces, Es, We),We}.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% wings_edge:select_region() does not work as I want with several loops
%% we know the faces though.
select_region(ELs, We) ->
    Es  = gb_sets:from_list([E || {Es,_} <- ELs, E <- Es]),
    Fs  = gb_sets:from_list([F || {_,Fs} <- ELs, F <- Fs]),
    ?dbg("Dissolve: ~w~n",[gb_sets:to_list(Fs)]),
    %% DBG = fun({Ess, Fss}) ->
    %% 		  FaceEs = wings_face:to_edges(Fss,We),
    %% 		  ?dbg("Face(~w) - Es ~w~n", [Fss,FaceEs -- Ess]),
    %% 		  ?dbg("Es - Face ~w~n", [Ess -- FaceEs])
    %% 	  end,
    %% [DBG(LI) || LI <- ELs],
    wings_edge:reachable_faces(Fs, Es, We).

weld(FsWe10, FsWe20) ->
    %% Store Fs1 and Fs2 in plugin so they are renumbered accordingly
    We0 = wings_we:merge(store_pst(FsWe10), store_pst(FsWe20)),
    %?dbg("After ~p: ~w~n",[We0#we.id,gb_trees:keys(We0#we.fs)]),
    {Faces, We1} = get_pst(We0),
    Weld = fun([F1,F2], WeAcc) -> do_weld(F1,F2,WeAcc) end,
    {#we{es=Etab} = We2, Es} = lists:foldl(Weld, {We1,[]}, Faces),
    Borders = ordsets:intersection(ordsets:from_list(Es),
                                   wings_util:array_keys(Etab)),
    %% Tesselate large border faces FIXME: only concave
    BorderFs = gb_sets:to_list(wings_face:from_edges(Borders, We2)),
    Fs = [Face || Face <- BorderFs, wings_face:vertices(Face, We2) > 5],
    We = wings_tesselation:quadrangulate(Fs, We2),
    {We, Borders}.

do_weld(Fa, Fb, {We0, Acc}) ->
    [Va|_] = wings_face:vertices_ccw(Fa, We0),
    Pos = wings_vertex:pos(Va, We0),
    Find = fun(Vb, _, _, Vs) ->
                   case wings_vertex:pos(Vb, We0) of
                       Pos -> [Vb|Vs];
                       _ -> Vs
                   end
           end,
    [Vb] = wings_face:fold(Find, [], Fb, We0),
    %% ?dbg("Fa ~p ~p (~p) Fb ~p ~p (~p)~n",
    %%      [Fa, Va, length(wings_face:vertices_ccw(Fa, We0)),
    %%       Fb, Vb, length(wings_face:vertices_ccw(Fb, We0))]),

    %% Bridge and collapse new edges
    We1 = wings_face_cmd:force_bridge(Fa, Va, Fb, Vb, We0),
    Es = wings_we:new_items_as_ordset(edge, We0, We1),
    We = lists:foldl(fun(E, W) -> wings_collapse:collapse_edge(E, W) end, We1, Es),
    %% Find selection
    BorderEdges = wings_face:to_edges([Fa,Fb], We0),
    {We, BorderEdges ++ Acc}.

order_loops(Fs, EsLoops0, We) ->
    %% ?dbg("Fs: ~w: ~w~n",[We#we.id, Fs]),
    %% [io:format("EL: ~w~n",[L]) || L <- EsLoops0],
    case length(Fs) =:= length(EsLoops0) of
        true -> ok;
        false -> throw(multiple_face_ints)
    end,
    FInfo = [{lists:sort(wings_face:to_edges([Face], We)), Face} || Face <- Fs],
    Loops = [lists:sort(Loop) || {Loop,_} <- EsLoops0],
    Pick = fun(Es) ->
                   {_, Face} = lists:keyfind(Es, 1, FInfo),
                   Face
           end,
    [Pick(Es) || Es <- Loops].

store_pst({Fs, #we{pst=Pst0}=We}) ->
    %?dbg("~p: ~w~n",[We#we.id, Fs]),
    {Data, _} = lists:mapfoldl(fun(F, I) -> {{I, F},I+1} end, 0, Fs),
    Pst = gb_trees:insert(?MODULE, gb_trees:from_orddict([{fs,Data}]), Pst0),
    We#we{pst=Pst}.

get_pst(#we{pst=Pst0}=We) ->
    Faces = lists:sort(gb_trees:get(?MODULE, Pst0)),
    Pst = gb_trees:delete(?MODULE, Pst0),
    {Faces, We#we{pst=Pst}}.

%% Pst callbacks
get_data(save, PData, Acc) ->
    FsInfo = gb_trees:get(fs,PData),
    {ok, [{plugin, {?MODULE, {fs,FsInfo}}}|Acc]};
get_data(_, _, Acc) -> Acc.

renumber(get_elem, {_Id, Face}) -> Face.
renumber(set_elem, {Id, _OldFace}, {value, NewFace}) ->
    %?dbg("ReN ~p: ~w => ~w~n", [Id, _OldFace, NewFace]),
    {Id, NewFace}.

merge_we([_,_]=Wes) ->
    try
        Lists = [FaceInfo || #we{pst=Pst} <- Wes, FaceInfo <- gb_trees:get(fs, gb_trees:get(?MODULE, Pst))],
        %?dbg("merge ~p~n",[Lists]),
        Rel = sofs:relation(Lists, [{id,data}]),
        Fam = sofs:relation_to_family(Rel),
        %% we ignore to make a gb_tree here..
        [Fs || {_, Fs} <- sofs:to_external(Fam)]
    catch _:Reason ->
            ?dbg("~p: ~p", [Reason, erlang:get_stacktrace()])
    end.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
make_verts(Loops, Vmap0, We0) ->
    {Vmap, We1} = cut_edges(lists:append(Loops), Vmap0, We0),
    Make = fun(Loop, We) ->
                   ?dbg("MakeEdgeLoop~n",[]), [io:format(" ~w~n", [L]) || L <- Loop],
                   make_edge_loop(Loop, Vmap, [], [], We)
           end,
    ELWe = {_,We} = lists:mapfoldl(Make, We1, Loops),
    ok = wings_we_util:validate(We),
    ELWe.

cut_edges(SE, Vmap, We0) ->
    WiEs = [{E,Vn} || #{op:=split_edge, e:=E, v:=Vn} <- SE],
    ECuts = sofs:to_external(sofs:relation_to_family(sofs:relation(WiEs, [{edge,vn}]))),
    lists:foldl(fun cut_edge/2, {Vmap, We0}, ECuts).

cut_edge({on_vertex, Vs}, {Vmap, #we{id=Id}=We}) ->
    {lists:foldl(fun(V, VM) ->
                         {Where, _Pos} = array:get(V, Vmap),
                         Vi = proplists:get_value(Id, Where),
                         array:set(V,Vi,VM)
                 end,
                 Vmap,Vs),
     We};
cut_edge({Edge, [V]}, {Vmap, #we{id=Id}=We0}) ->
    {Where,Pos} = array:get(V, Vmap),
    case proplists:get_value(Id, Where) of
        undefined ->
            {We, NewV} = wings_edge:fast_cut(Edge, Pos, We0),
            {array:set(V, NewV, Vmap), We};
        Vi ->
            {array:set(V, Vi, Vmap), We0}
    end;
cut_edge({Edge, Vs}, {Vmap0, #we{es=Etab}=We0}) ->
    #edge{vs=VS} = array:get(Edge, Etab),
    P1 = wings_vertex:pos(VS,We0),
    C = fun(V) ->
                Pos = vmap_pos(V, Vmap0),
                Dist2 = e3d_vec:dist_sqr(Pos, P1),
                {Dist2, V, Pos}
        end,
    VsPos = lists:sort([C(V) || V <- Vs]),
    {We,_,Vmap} = lists:foldl(fun({_,V,Pos}, {WE, E, Vm}) ->
                                      {We, New} = wings_edge:fast_cut(E, Pos, WE),
                                      {We, New, array:set(V, New, Vm)}
                              end, {We0, Edge, Vmap0}, VsPos),
    {Vmap, We}.

make_edge_loop([#{op:=split_edge}=F|_]=Loop, Vmap, EL, IFs, We) ->
    make_edge_loop_1(Loop, F, Vmap, EL, IFs, We);
make_edge_loop(Loop, Vmap, EL, IFs, We) ->
    %% Start with split_edge
    case lists:splitwith(fun(#{op:=Op}) -> Op =:= split_face end, Loop) of
        {FSs, []} -> split_face(FSs, Vmap, EL, We);
        {FSs, Edges} -> make_edge_loop(Edges++FSs, Vmap, EL, IFs, We)
    end.

make_edge_loop_1([V1], V1, _, EL, IFs, We) ->
    {{EL, IFs}, We};
make_edge_loop_1([#{op:=split_edge}=V1],#{op:=split_edge}=V2, Vmap, EL, IFs, We0) ->
    {{We, New}, Face} = connect_verts(V1,V2,Vmap, We0),
    {{[New|EL], Face ++ IFs}, We};
make_edge_loop_1([#{op:=split_edge}=V1|[#{op:=split_edge}=V2|_]=Rest], Last, Vmap, EL, IFs, We0) ->
    {{We, New}, Face} = connect_verts(V1,V2,Vmap, We0),
    make_edge_loop_1(Rest, Last, Vmap, [New|EL], Face ++ IFs, We);
make_edge_loop_1([#{op:=split_edge}=V1|Splits], Last, Vmap, EL0, IFs, We0) ->
    {FSs, [V2|_]=Rest} =
        case lists:splitwith(fun(#{op:=Op}) -> Op =:= split_face end, Splits) of
            {FS, []} -> {FS, [Last]};
            {FS, Rs} -> {FS, Rs}
        end,
    case edge_exists(V1,V2,Vmap,We0) of
        [] -> %% Standard case
            ?dbg("Connect: ~w[~w] ~w[~w]~n",[maps:get(v,V1), array:get(maps:get(v,V1),Vmap),
                                             maps:get(v,V2), array:get(maps:get(v,V2), Vmap)]),
            {{We1, Edge}, Face} = connect_verts(V1,V2,{EL0,FSs},Vmap,We0),
            ok = wings_we_util:validate(We1),
            ?dbg("new edge ~w face ~w~n",[Edge, Face]),
            {EL1,Vmap1,We2} = make_face_vs(FSs, V1, Edge, Vmap, We1),
            make_edge_loop_1(Rest, Last, Vmap1, EL1++EL0, Face++IFs, We2);
        [{Edge,_F1,_F2}] ->
            ?dbg("Ignore ~w ~w edge ~p~n",[maps:get(v,V1),maps:get(v,V2),Edge]),
            {EL1,Vmap1,We} = make_face_vs(FSs, V1, Edge, Vmap, We0),
            make_edge_loop_1(Rest, Last, Vmap1, EL1++EL0, IFs, We)
    end.

edge_exists(#{v:=V10},#{v:=V20},Vmap,We) ->
    wings_vertex:edge_through(array:get(V10,Vmap),array:get(V20,Vmap),We).

connect_verts(V1, V2, Vmap, We) ->
    connect_verts(V1, V2, [], Vmap, We).
connect_verts(V1, V2, RefPoints, Vmap, #we{vp=Vtab}=We) ->
    {WeV1,WeV2,Face,OtherN} = pick_face(RefPoints, V1,V2, Vmap, We),
    case wings_vertex:edge_through(WeV1,WeV2,Face,We) of
        none ->
            ?dbg("~w: ~w ~w in ~w ~s~n",[We#we.id, WeV1, WeV2, Face, e3d_vec:format(OtherN)]),
            N = wings_face:normal(Face, We),
            Dir = e3d_vec:cross(N,e3d_vec:sub(array:get(WeV1,Vtab),array:get(WeV2,Vtab))),
	    ?dbg("Swap: ~p~n", [0 >= e3d_vec:dot(OtherN, Dir)]),
            case 0 >= e3d_vec:dot(OtherN, Dir) of
                true  -> {wings_vertex:force_connect(WeV1,WeV2,Face,We), [Face]};
                false -> {wings_vertex:force_connect(WeV2,WeV1,Face,We), [Face]}
            end;
        Edge when RefPoints =:= [] ->
            ?dbg("Skip ~p ~p~n",[Edge,Face]),
            {{We, Edge}, []}
    end.

pick_face([], #{v:=V1,o_n:=N1}, #{v:=V2,o_n:=N2}, Vmap, We) ->
    WeV1 = array:get(V1, Vmap),
    WeV2 = array:get(V2, Vmap),
    true = is_integer(WeV1), true = is_integer(WeV2), %% Assert
    OtherN = case e3d_vec:norm(e3d_vec:average(N1,N2)) of
                 {0.0,0.0,0.0} -> error({N1,N2});
                 ON -> ON
             end,
    case [Face || {Face, [_,_]} <- wings_vertex:per_face([WeV1,WeV2],We)] of
        [Face] ->
            {WeV1,WeV2,Face,OtherN};
        [Face|_] = _Fs ->
            {WeV1,WeV2,Face,OtherN}
    end;
pick_face({Edges,[#{v:=V3, o_n:=N}|_]=Refs}, #{v:=V1,fs:={LF0,RF0}}, #{v:=V2},
          Vmap, #we{es=Etab, vp=Vtab}=We) ->
    WeV1 = array:get(V1, Vmap),
    WeV2 = array:get(V2, Vmap),
    All = wings_vertex:per_face([WeV1,WeV2],We),
    case [Face || {Face, [_,_]} <- All] of
        [Face] ->
            {WeV1,WeV2,Face, N};
        [F1,F2] ->
	    V3P  = vmap_pos(V3, Vmap),
	    Dir = e3d_vec:sub(V3P,array:get(WeV1,Vtab)),
            [Edge|_] = Edges,
            ?dbg("id:~p V=~p Fs: ~p ~p~n", [We#we.id, WeV1, F1,F2]),
            {LF,RF,_} = _Dbg =
                case array:get(Edge, Etab) of
                    #edge{vs=WeV1, lf=F1, rf=F2} -> {F1,F2,?LINE};
                    #edge{vs=WeV1, lf=F2, rf=F1} -> {F2,F1,?LINE};
                    #edge{ve=WeV1, lf=F1, rf=F2} -> {F1,F2,?LINE};
                    #edge{ve=WeV1, lf=F2, rf=F1} -> {F2,F1,?LINE}
                end,
            Face = case pick_ref_face(Refs, undefined) of
                       LF0 -> LF;
                       RF0 -> RF
                   end,
            ?dbg("pick ~p in ~w ~w => ~p~n", [pick_ref_face(Refs, undefined), [LF0,RF0], _Dbg, Face]),
            {WeV1,WeV2,Face,Dir}
    end.

pick_ref_face([#{f:=F}|Ss], undefined) ->
    pick_ref_face(Ss, F);
pick_ref_face([#{f:=F}|Ss], F) ->
    pick_ref_face(Ss, F);
pick_ref_face([], F) -> F.

make_face_vs([_]=Ss, _Vs, Edge, Vmap, We) ->
    make_face_vs_1(Ss, Edge, Vmap, [Edge], We);
make_face_vs(Ss, #{v:=Vs0}, Edge, Vmap, #we{es=Etab}=We) ->
    Vs = array:get(Vs0, Vmap),
    case array:get(Edge, Etab) of
        #edge{vs=Vs} -> make_face_vs_1(Ss, Edge, Vmap, [Edge], We);
        #edge{ve=Vs} ->
            {EL1,VM1,WE1} = make_face_vs_1(lists:reverse(Ss), Edge, Vmap, [Edge], We),
            {lists:reverse(EL1),VM1,WE1}
    end.

make_face_vs_1([#{op:=split_face,v:=V}|Ss], Edge, Vmap, EL, We0) ->
    Pos = vmap_pos(V, Vmap),
    {We, New} = wings_edge:fast_cut(Edge, Pos, We0),
    make_face_vs_1(Ss, New, array:set(V, New, Vmap), [New|EL], We);
make_face_vs_1([], _, Vmap, EL, We) ->
    {EL, Vmap, We}.

split_face([#{v:=V1},#{v:=V2}|_]=Fs, Vmap, EL, We0) ->
    Face = pick_ref_face(Fs, undefined),
    NumberOfNew = length(Fs),
    true = NumberOfNew > 2, %% Otherwise something is wrong
    We1 = wings_extrude_face:faces([Face], We0),
    FVs = wings_face:vertices_ccw(Face, We1),
    [P1,P2|_] = FPos = wings_face:vertex_positions(Face, We1),
    Zipped = lists:zip(FVs, FPos),
    NumberOfOld = length(FVs),
    if
	NumberOfOld =:= NumberOfNew ->
	    %?dbg("~p ~n",[Zipped]),
	    KD3 = e3d_kd3:from_list(Zipped),
	    Center = e3d_vec:average(FPos),
	    %?dbg("~p ~p~n", [Center, FPos]),
	    D1 = e3d_vec:normal(P1,P2,Center),
	    P3 = vmap_pos(V1, Vmap),
	    P4 = vmap_pos(V2, Vmap),
	    D2 = e3d_vec:normal(P3,P4,Center),
	    %io:format("~p ~p => ~p~n", [D1,D2,e3d_vec:dot(D1,D2)]),
	    Ordered = case e3d_vec:dot(D1,D2) > 0 of
			  true -> FVs;
			  false -> lists:reverse(FVs)
		      end,
	    {{First,_}, _} = e3d_kd3:take_nearest(P3, KD3),
	    %?dbg("Search ~p~n", [P3]),
	    {VL1,VL2} = lists:splitwith(fun(V) when V =:= First -> false; (_) -> true end,
					Ordered),
	    PosL = [vmap_pos(Vi, Vmap)|| #{v:=Vi} <- Fs],
	    %?dbg("~w ~w ~p~n",[First,VL2++VL1,PosL]),
	    Vs = lists:zip(VL2++VL1, PosL),
	    Vtab = lists:foldl(fun({V,Pos}, Vtab) -> array:set(V, Pos, Vtab) end,
			       We1#we.vp, Vs),
            cleanup_edges(FVs, [V||{V,_}<-Vs], Face, EL, We1#we{vp=Vtab});
	NumberOfOld > NumberOfNew ->
	    KD3 = e3d_kd3:from_list(Zipped),
	    {Vs,_} = lists:mapfoldl(fun(#{v:=Vi}, Tree0) ->
					    Pos = vmap_pos(Vi, Vmap),
					    {{V,_}, Tree} = e3d_kd3:take_nearest(Pos, Tree0),
					    {{V,Pos},Tree}
				    end, KD3, Fs),
	    Vtab = lists:foldl(fun({V,Pos}, Vtab) -> array:set(V, Pos, Vtab) end,
			       We1#we.vp, Vs),
            cleanup_edges(FVs, [V||{V,_}<-Vs], Face, EL, We1#we{vp=Vtab});
	true ->
	    KD3 = e3d_kd3:from_list([{FS, vmap_pos(Vi, Vmap)} || #{v:=Vi}=FS <- Fs]),
	    {Vs,_} = lists:mapfoldl(fun({V, Old}, Tree0) ->
					    {{FS,Pos}, Tree} = e3d_kd3:take_nearest(Old, Tree0),
					    {{V,Pos,FS},Tree}
				    end, KD3, Zipped),
	    Vtab = lists:foldl(fun({V, Pos, _}, Vtab) -> array:set(V, Pos, Vtab) end,
			       We1#we.vp, Vs),
	    Fs1 = lists:map(fun(FS) -> case lists:keyfind(FS, 3, Vs) of
					   false -> FS;
					   {_,_,_} -> FS#{op:=split_edge, o_n=>ignore}
				       end
			    end, Fs),
	    Vmap1 = lists:foldl(fun({V, _, #{v:=Vi}}, Map) -> array:set(Vi, V, Map) end,
				Vmap, Vs),
	    {{Es,_},We} = make_edge_loop(Fs1, Vmap1, EL, [], We1#we{vp=Vtab}),
            {{Es,[Face]},We}
    end.

cleanup_edges(FVs, Used, Face, EL, We) ->
    %% Start with a used vertex
    {Vs1,Vs0} = lists:splitwith(fun(V) -> not lists:member(V, Used) end, FVs),
    %% ?dbg("~p ~p~n",[EL, Vs0++Vs1]),
    cleanup_edges(Vs0++Vs1, false, hd(Vs0), [], Used, Face, EL, We).

cleanup_edges([V1|[V2|Vs]=Vs0], Connect, Last, Drop, Used, Face, EL, We0) ->
    case lists:member(V2, Used) of
        true when Connect ->
            {We, New} = wings_vertex:force_connect(V2,V1,Face,We0),
            cleanup_edges(Vs0, false, Last, Drop, Used, Face, [New|EL], We);
        true ->
            Edge = wings_vertex:edge_through(V1,V2,Face,We0),
            cleanup_edges(Vs0, false, Last, Drop, Used, Face, [Edge|EL], We0);
        false ->
            cleanup_edges([V1|Vs], true, Last, [V2|Drop], Used, Face, EL, We0)
    end;
cleanup_edges([V1], Connect, Last, Drop, _Used, Face, EL0, We0) ->
    {EL,We2} = case Connect of
                   true ->
                       {We1, Edge} = wings_vertex:force_connect(Last,V1,Face,We0),
                       {[Edge|EL0],We1};
                   false ->
                       Edge = wings_vertex:edge_through(V1,Last,Face,We0),
                       {[Edge|EL0],We0}
               end,
    % ?dbg("drop vs ~p~n",[Drop]),
    Es = wings_edge:from_vs(Drop, We2),
    %% ?dbg("drop es ~p ~w~n",[Es, EL]),
    We3 = wings_edge:dissolve_edges(Es, We2),
    ok = wings_we_util:validate(We3),
    {{EL, [Face]},We3}.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filter_tri_edges({L1,L2}, We1,We2) ->
    Loop = lists:zip(L1,L2),
    %% ?dbg("~p~n",[?FUNCTION_NAME]), [io:format("1: ~w~n2: ~w~n~n",[K1,K2])||{K1,K2}<-Loop],
    Res = filter_tri_edges_1(Loop,We1,We2),
    %% ?dbg("after ~p~n",[?FUNCTION_NAME]), [io:format("1: ~w~n2: ~w~n~n",[K1,K2])||{K1,K2}<-Res],
    lists:unzip(Res).

filter_tri_edges_1([{#{v:=V}=V1,U1}, {#{v:=V}=V2, U2}|Vs],We1,We2) ->
    %% Remove edges to it self (loops)
    filter_tri_edges_1([{filter_edge(V1,V2),filter_edge(U1,U2)}|Vs],We1,We2);

filter_tri_edges_1([{#{op:=split_edge,e:=none},#{op:=split_edge, e:=none}}|Vs],We1,We2) ->
    filter_tri_edges_1(Vs,We1,We2);
filter_tri_edges_1([{#{op:=split_edge,e:=none,f:=F}=V1,#{op:=Op}=V2}|Vs],We1,We2) ->
    case Op of
        split_face -> filter_tri_edges_1(Vs,We1,We2);
        split_edge ->
            case skip_tess_edge(wings_face:normal(F,We1), V2, We2) of
                true -> filter_tri_edges_1(Vs,We1,We2);
                false -> [{edge_to_face(V1), V2}|filter_tri_edges_1(Vs,We1,We2)]
            end
    end;
filter_tri_edges_1([{#{op:=Op}=V1,#{op:=split_edge,e:=none,f:=F}=V2}|Vs],We1,We2) ->
    case Op of
        split_face -> filter_tri_edges_1(Vs,We1,We2);
        split_edge ->
            case skip_tess_edge(wings_face:normal(F,We2), V1, We1) of
                true -> filter_tri_edges_1(Vs,We1,We2);
                false -> [{V1,edge_to_face(V2)}|filter_tri_edges_1(Vs,We1,We2)]
            end
    end;
filter_tri_edges_1([V|Vs],We1,We2) ->
    [V|filter_tri_edges_1(Vs,We1,We2)];
filter_tri_edges_1([],_We1,_We2) -> [].

filter_edge(_, #{op:=split_edge, e:=Edge}=V2) when Edge =/= none -> V2;
filter_edge(V1,_) -> V1.

skip_tess_edge(N, #{e:=Edge}=EC, #we{es=Etab}=We) ->
    #edge{vs=VS,ve=VE} = array:get(Edge,Etab),
    Dir = e3d_vec:sub(wings_vertex:pos(VS, We),wings_vertex:pos(VE,We)),
    ?dbg("~p~n",[EC]),
    ?dbg("~p ~p => ~p~n", [N, Dir, e3d_vec:dot(N, Dir)]),
    abs(e3d_vec:dot(N, Dir)) < 0.1.

edge_to_face(#{op:=split_edge}=Orig) ->
    Orig#{op=>split_face}.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% We need to build the cycle our selfves since the edges may not be directed
%% in the correct direction. Also the V is new Vs on edges and there maybe
%% several new V on the same wings edge.

build_vtx_loops(Edges, _Acc) ->
    G = make_lookup_table(Edges),
    Comps = digraph_utils:components(G),
    %% ?dbg("Cs: ~w~n",[Comps]),
    Res = [build_vtx_loop(C, G) || C <- Comps],
    [] = digraph:edges(G), %% Assert that we have completed all edges
    digraph:delete(G),
    Res.

make_lookup_table(Edges) ->
    G = digraph:new(),
    Add = fun(#{p1:={_,V},p2:={_,V}}) ->
                  ignore; %% Loops
             (#{p1:={_,P1},p2:={_,P2}}=EI) ->
                  digraph:add_vertex(G, P1),
                  digraph:add_vertex(G, P2),
                  case edge_exists(G,P1,P2) of
                      false -> digraph:add_edge(G, P1, P2, EI);
                      true -> ok
                  end
          end,
    _ = [Add(EI) || EI <- Edges],
    G.

build_vtx_loop([V|_Vs], G) ->
    case build_vtx_loop(V, G, []) of
        {V, Acc} -> Acc;
        {_V, _Acc} ->
            ?dbg("V ~p => ~p~n",[V,_Acc]),
            ?dbg("Last ~p~n",[_V]),
            error(incomplete_edge_loop)
    end.

build_vtx_loop(V0, G, Acc) ->
    case [digraph:edge(G, E) || E <- digraph:edges(G, V0)] of
        [] -> {V0, Acc};
        Es ->
            {Edge, Next, Ei} = pick_edge(Es, V0, undefined),
            %?dbg("~p in ~P~n => ~p ~n",[V0, Es, 10, Next]),
            digraph:del_edge(G, Edge),
            build_vtx_loop(Next, G, [Ei,V0|Acc])
    end.

edge_exists(G,V1,V2) ->
    lists:member(V2, digraph:out_neighbours(G, V1)) orelse
        lists:member(V1, digraph:out_neighbours(G, V2)).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pick_edge([{E,V,V,Ei}|_], V, _Best) ->
    {E, V, Ei}; %% Self cyclic pick first
pick_edge([{E,V,N,Ei}|R], V, _Best) ->
    pick_edge(R, V, {E,N,Ei});
pick_edge([{E,N,V,Ei}|R], V, _Best) ->
    pick_edge(R, V, {E,N,Ei});
pick_edge([], _, Best) -> Best.

split_loop([Last|Loop], Vmap, We) ->
    split_loop(Loop, Last, Vmap, We, []).

split_loop([V1,E|Loop], Last, Vmap, We, Acc) when is_integer(V1) ->
%    ?dbg("~p: ~p in ~p~n",[(element(1,We))#we.id,V1,E]),
    Vertex = vertex_info(E, V1, Vmap, We),
    split_loop(Loop, Last, Vmap, We, [Vertex|Acc]);
split_loop([V1], E, Vmap, We, Acc) ->
%    ?dbg("~p: ~p in ~p~n",[(element(1,We))#we.id,V1,E]),
    Vertex = vertex_info(E, V1, Vmap, We),
    lists:reverse([Vertex|Acc]).

vertex_info(#{mf1:={O1,F1}, mf2:={O2,F2}, other:={O3,F3},
              p1:={{_, {A1,B1}=Edge1},V0},
              p2:={{_, {A2,B2}=Edge2},V0}},
            V0, Vmap, {#we{id=Id}=We,OWe}) ->
    if O1 =:= Id ->
            ON = if O2=:=Id -> wings_face:normal(F3, OWe);
                    true -> wings_face:normal(F2, OWe)
                end,
            Edge = wings_vertex:edge_through(A1,B1,F1,We),
            Fs = edge_faces(Edge,F1,We),
            #{op=>split_edge, o=>Id, f=>F1, e=>Edge, v=>V0, vs=>Edge1, fs=>Fs, o_n=>ON};
       O2 =:= Id ->
            ON = if O1=:=Id -> wings_face:normal(F3, OWe);
                    true -> wings_face:normal(F1, OWe)
                 end,
            Edge = wings_vertex:edge_through(A2,B2,F2,We),
            Fs = edge_faces(Edge,F1,We),
            #{op=>split_edge, o=>Id, f=>F2, e=>Edge, v=>V0, vs=>Edge2, fs=>Fs, o_n=>ON};
       O3 =:= Id ->
            ON = wings_face:normal(F1, OWe),
            check_if_edge(#{op=>split_face, o=>Id, f=>F3, v=>V0, o_n=>ON}, Vmap, We)
    end;
vertex_info(#{mf1:={O1,F1}, mf2:={O2,F2}, other:={O3,F3}, p1:={{_, {A,B}=Edge0},V0}}, V0,
            Vmap, {#we{id=Id}=We,OWe}) ->
    if O1 =:= Id ->
            ON = if O2=:=Id -> wings_face:normal(F3, OWe);
                   true -> wings_face:normal(F2, OWe)
                end,
            Edge = wings_vertex:edge_through(A,B,F1,We),
            Fs = edge_faces(Edge,F1,We),
            #{op=>split_edge, o=>Id, f=>F1, e=>Edge, v=>V0, vs=>Edge0, fs=>Fs, o_n=>ON};
       O2 =:= Id ->
            ON = wings_face:normal(F1, OWe),
            check_if_edge(#{op=>split_face, o=>Id, f=>F2, v=>V0, o_n=>ON}, Vmap, We);
       O3 =:= Id ->
            ON = wings_face:normal(F1, OWe),
            check_if_edge(#{op=>split_face, o=>Id, f=>F3, v=>V0, o_n=>ON}, Vmap, We)
    end;
vertex_info(#{mf2:={O1,F1}, mf1:={O2,F2}, other:={O3,F3}, p2:={{_, {A,B}=Edge0},V0}}, V0,
            Vmap, {#we{id=Id}=We,OWe}) ->
    if O1 =:= Id ->
            ON = if O2=:=Id -> wings_face:normal(F3, OWe);
                    true -> wings_face:normal(F2, OWe)
                 end,
            Edge = wings_vertex:edge_through(A,B,F1,We),
            Fs = edge_faces(Edge,F1,We),
            #{op=>split_edge, o=>Id, f=>F1, e=>Edge, v=>V0, vs=>Edge0, fs=>Fs, o_n=>ON};
       O2 =:= Id ->
            ON = wings_face:normal(F1, OWe),
            check_if_edge(#{op=>split_face, o=>Id, f=>F2, v=>V0, o_n=>ON}, Vmap, We);
       O3 =:= Id ->
            ON = wings_face:normal(F1, OWe),
            check_if_edge(#{op=>split_face, o=>Id, f=>F3, v=>V0, o_n=>ON}, Vmap, We)
    end.

edge_faces(none,F1, _We) ->
    {F1,F1};
edge_faces(Edge,_F1, #we{es=Etab}) ->
    #edge{lf=LF,rf=RF} = array:get(Edge, Etab),
    {LF,RF}.

check_if_edge(#{f:=F, v:=V}=SF, Vmap, #we{id=Id, vp=Vtab, es=Etab}=We) ->
    {Where, Pos} = array:get(V, Vmap),
    Find = fun(_,Edge,#edge{vs=V1,ve=V2},Acc) ->
                   V1P = array:get(V1, Vtab),
                   V2P = array:get(V2, Vtab),
                   case e3d_vec:line_dist_sqr(Pos, V1P, V2P) < ?EPSILON of
                       true -> [{{V1,V2}, Edge}|Acc];
                       false -> Acc
                   end
           end,
    Es = wings_face:fold(Find, [], F, We),
    case {proplists:get_value(Id, Where), Es} of
        {undefined, []} -> SF;
        {undefined, [{Vs,Edge}]} ->
            #edge{lf=LF,rf=RF} = array:get(Edge, Etab),
            SF#{op:=split_edge, e=>Edge, vs=>Vs, fs=>{LF,RF}};
        {WeV, [{_,Edge}|_]} ->
            #edge{lf=LF,rf=RF} = array:get(Edge, Etab),
            SF#{op:=split_edge, e=>on_vertex, fs=>{LF,RF}, vs=>WeV}
    end.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

make_bvh(_, #we{id=Id, fs=Fs0}=We0, Bvhs) ->
    Fs = gb_trees:keys(Fs0),
    {Vtab,Ts} = triangles(Fs, We0),
    Get = fun({verts, Face}) -> element(1, array:get(Face, Ts));
	     (verts) -> Vtab;
	     (meshId) -> Id
	  end,
    Bvh = e3d_bvh:init([{array:size(Ts), Get}]),
    [#{id=>Id,map=>Ts,bvh=>Bvh,fs=>[],we=>We0}|Bvhs].

make_vmap(ReEI0, #we{id=Id1, vp=Vtab1}, #we{id=Id2, vp=Vtab2}) ->
    {I0,L1} = lists:foldl(fun({N, Pos}, {I, Acc}) ->
                                  {I+1, [{{I, [{Id1,N}]}, Pos}|Acc]}
                          end, {0, []}, array:sparse_to_orddict(Vtab1)),
    Tree0 = e3d_kd3:from_list(L1),
    {I1,Tree1} = add_vtab(Vtab2, I0, Id2, Tree0),
    make_vmap(ReEI0, Tree1, I1, []).

add_vtab(Vtab, I0, Id, Tree) ->
    Add = fun(N, Pos, {I,Acc}) ->
                  {{IF,V1},P1} = Obj = e3d_kd3:nearest(Pos, Acc),
                  New = {Id,N},
                  case e3d_vec:dist_sqr(Pos, P1) < ?EPSILON of
                      true  -> {I, e3d_kd3:update(Obj, {IF,[New|V1]}, Acc)};
                      false -> {I+1, e3d_kd3:enter(Pos, {I, [New]}, Acc)}
                  end
          end,
    array:sparse_foldl(Add, {I0,Tree}, Vtab).

make_vmap([#{p1:=P10, p2:=P20}=E|R], T0, N0, Acc) ->
    {P1, N1, T1} = vmap(P10, N0, T0),
    {P2, N2, T2} = vmap(P20, N1, T1),
    make_vmap(R, T2, N2, [E#{p1:=P1,p2:=P2}|Acc]);
make_vmap([], T, _, Acc) ->
    OrdD = [{N,{Where,Pos}} || {{N, Where}, Pos} <- lists:sort(e3d_kd3:to_list(T))],
    {array:from_orddict(OrdD), Acc}.

vmap({Where, Pos}, N, Tree) ->
    {{I, _V1}, P1} = e3d_kd3:nearest(Pos, Tree),
    case e3d_vec:dist_sqr(Pos, P1) < ?EPSILON of
        true  -> {{Where, I}, N, Tree};
        false -> {{Where, N}, N+1, e3d_kd3:enter(Pos, {N, []}, Tree)}
    end.

vmap_pos(N, Vmap) ->
    {_Where, Pos} = array:get(N, Vmap),
    Pos.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

remove(#{id:=Id}, List) ->
    Map = mapsfind(Id, id, List),
    lists:delete(Map, List).

remove(H1, H2, List) ->
    remove(H1, remove(H2, List)).

remap(#{mf1:=MF10,mf2:=MF20,p1:={Pos1,E11,E12}, p2:= {Pos2, E21, E22}, other:=Other},
      #{id:=Id1,map:=M1}, #{id:=Id2,map:=M2}) ->
    MF1 = remap_1(MF10, Id1, M1, Id2, M2),
    MF2 = remap_1(MF20, Id1, M1, Id2, M2),
    Oth = remap_1(Other, Id1, M1, Id2, M2),
    EId1 = {element(1,MF1), order(E11,E12)},
    EId2 = {element(1,MF2), order(E21,E22)},
    case Id1 < Id2 of
        true  -> #{mf1=>MF1, mf2=>MF2, p1=>{EId1, Pos1}, p2=>{EId2, Pos2}, other=>Oth};
        false -> #{mf1=>MF2, mf2=>MF1, p1=>{EId2, Pos2}, p2=>{EId1, Pos1}, other=>Oth}
    end.

remap_1({Id, TriFace}, Id, M1, _Id2, _M2) ->
    {_, Face} = array:get(TriFace, M1),
    {Id, Face};
remap_1({Id, TriFace}, _Id, _M1, Id, M2) ->
    {_, Face} = array:get(TriFace, M2),
    {Id, Face}.

order(V1, V2) when V1 < V2 ->  {V1,V2};
order(V2, V1) -> {V1,V2}.

%% Use wings_util:mapsfind
mapsfind(Value, Key, [H|T]) ->
    case H of
	#{Key:=Value} -> H;
	_ -> mapsfind(Value, Key, T)
    end;
mapsfind(_, _, []) -> false.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

triangles(Fs, We0) ->
    {#we{vp=Vtab}, Ts} = lists:foldl(fun triangle/2, {We0,[]}, Fs),
    {Vtab, array:from_list(Ts)}.

triangle(Face, {We, Acc}) ->
    Vs = wings_face:vertices_ccw(Face, We),
    case length(Vs) of
	3 -> {We, [{list_to_tuple(Vs), Face}|Acc]};
	4 -> tri_quad(Vs, We, Face, Acc);
	_ -> tri_poly(Vs, We, Face, Acc)
    end.

tri_quad([Ai,Bi,Ci,Di] = Vs, #we{vp=Vtab}=We, Face, Acc) ->
    [A,B,C,D] = VsPos = [array:get(V, Vtab) || V <- Vs],
    N = e3d_vec:normal(VsPos),
    case wings_tesselation:is_good_triangulation(N, A, B, C, D) of
	true  -> {We, [{{Ai,Bi,Ci}, Face}, {{Ai,Ci,Di}, Face}|Acc]};
	false -> {We, [{{Ai,Bi,Di}, Face}, {{Bi,Ci,Di}, Face}|Acc]}
    end.

tri_poly(Vs, #we{vp=Vtab}=We, Face, Acc0) ->
    VsPos = [array:get(V, Vtab) || V <- Vs],
    N = e3d_vec:normal(VsPos),
    {Fs0, Ps0} = wings_gl:triangulate(N, VsPos),
    ?dbg("~p ~w~n",[Face, Vs]),
    io:format("~p~n", [Fs0]),
    [io:format("~s ",[e3d_vec:format(V)]) || V <- VsPos],
    io:format("~n", []),
    [io:format("~s ",[e3d_vec:format(V)]) || V <- Ps0],
    io:format("~n", []),
    Index = array:size(Vtab),
    {TessVtab, Is} = renumber_and_add_vs(Ps0, Vs, Index, Vtab, []),
    io:format("~p~n", [Is]),
    Fs = lists:foldl(fun({A,B,C}=_F, Acc) ->
			     F = {element(A,Is), element(B, Is), element(C,Is)},
                             io:format("~p => ~p~n", [_F, F]),
			     [{F, Face}|Acc]
		     end, Acc0, Fs0),
    {We#we{vp=TessVtab}, Fs}.

renumber_and_add_vs([_|Ps], [V|Vs], Index, Vtab, Acc) ->
    renumber_and_add_vs(Ps, Vs, Index, Vtab, [V|Acc]);
renumber_and_add_vs([Pos|Ps], [], Index, Vtab0, Acc) ->
    Vtab = array:set(Index, Pos, Vtab0),
    renumber_and_add_vs(Ps, [], Index+1, Vtab, [Index|Acc]);
renumber_and_add_vs([], [], _, Vtab, Vs) ->
    {Vtab, list_to_tuple(lists:reverse(Vs))}.

