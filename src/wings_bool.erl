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
    IntScts = find_intersects(MapBvh, MapBvh, []),
    MFs = lists:append([MFL || #{es:=MFL} <- IntScts]),
    Sel0 = sofs:relation(MFs, [{id,data}]),
    Sel1 = sofs:relation_to_family(Sel0),
    Sel2 = sofs:to_external(Sel1),
    Sel = [{Id,gb_sets:from_list(L)} || {Id,L} <- Sel2],
    %% ?dbg("~p~n", [Sel2]),
    Upd = fun(#{we:=Wes}, Sh) ->
                  case Wes of
                      We=#we{id=Id} ->
                          gb_trees:update(Id, We, Sh);
                      {We1=#we{id=Id1},We2=#we{id=Id2, vp=Vtab0}} ->
                          Sh2 = gb_trees:update(Id1, We1, Sh),
                          Vtab = array:sparse_map(fun(_, V) -> e3d_vec:add({0.1,0.1,0.0},V) end, Vtab0),
                          gb_trees:update(Id2, We2#we{vp=Vtab}, Sh2)
                  end
          end,
    Sh = lists:foldl(Upd, Sh0, IntScts),
    wings_sel:set(edge, Sel, St0#st{shapes=Sh}).

find_intersects([Head|Tail], [_|_] = Alts, Acc) ->
    case find_intersect(Head, Alts) of
	false ->
	    find_intersects(Tail,Alts, Acc);
	{H1,Merge} ->
	    NewTail = [Merge|remove(H1, Head, Tail)],
	    NewAlts = [Merge|remove(H1, Head, Alts)],
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
    Loops = [filter_tri_edges(Loop) || Loop <- lists:zip(L10,L20)],
    {L1,L2} = lists:unzip(Loops),
    {Es1, We1N0} = make_verts(L1, Vmap, We1),
    {Es2, We2N0} = make_verts(L2, Vmap, We2),
    GetFaces = fun(Id, Ls) -> [{Id,F} || Loop <- Ls, #{f:=F} <- Loop] end,
    MFs = GetFaces(Id1, L1) ++ GetFaces(Id2, L2),
    GetEdges = fun(Id, Es) -> [{Id,E} || Loop <- Es, E <- Loop] end,
    MEs = GetEdges(Id1, Es1) ++ GetEdges(Id2, Es2),

    DRes1 = {_,We1N} = dissolve_faces_in_edgeloops(Es1, We1N0),
    DRes2 = {_,We2N} = dissolve_faces_in_edgeloops(Es2, We2N0),

    weld(DRes1, DRes2),

    I1#{id:=min(Id1,Id2), es=>MEs, fs:=MFs, we:={We1N,We2N}}.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dissolve_faces_in_edgeloops(Es, We0) ->
    [{_, Fs}] = wings_edge:select_region(lists:append(Es), We0),
    %% verify or invert faces
    We = wings_dissolve:faces(Fs, We0),
    Faces0 = wings_we:new_items_as_ordset(face, We0, We),
    Faces = Faces0, %% reorder with Es so they are in Es order
    {Faces, We}.

weld({Fs1,We10}, {Fs2, We20}) ->
    %% Store Fs1 and Fs2 in plugin so they are renumbered accordingly
    We0 = wings_we:merge(We10, We20),
    %% Fetch Fs1 and Fs2 from pst
    lists:foldl(fun({F1,F2}, We) -> do_weld(F1,F2,We) end, We0, lists:zip(Fs1,Fs2)).

do_weld(_F1, _F2, We) ->
    % Pick Va and Vb see wings_body:try_weld_1
    %% wings_face_cmd:force_bridge(Fa, Va, Fb, Vb, We0),
    We.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
make_verts(Loops, Vmap0, We0) ->
    All = lists:append(Loops),
    SE = [SP || #{op:=split_edge}=SP <- All],
    {Vmap, We1} = cut_edges(SE, Vmap0, We0),
    Make = fun(Loop, We) -> make_edge_loop(Loop, Vmap, [], We) end,
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
    #edge{vs=VS,ve=VE} = array:get(Edge, Etab),
    Pt1 = wings_vertex:pos(VS,We0),
    Pt2 = wings_vertex:pos(VE,We0),
    [{_,Pt3},{_,Pt4}|_] = VsPos0 = [{V,vmap_pos(V, Vmap0)} || V <- Vs],
    VsPos = case e3d_vec:dot(e3d_vec:sub(Pt2, Pt1), e3d_vec:sub(Pt4, Pt3)) > 0 of
                true -> VsPos0;
                false -> lists:reverse(VsPos0)
            end,
    {We,_,Vmap} = lists:foldl(fun({V,Pos}, {WE, E, Vm}) ->
                                      {We, New} = wings_edge:fast_cut(E, Pos, WE),
                                      {We, New, array:set(V, New, Vm)}
                              end, {We0, Edge, Vmap0}, VsPos),
    {Vmap, We}.

make_edge_loop([#{op:=split_edge}=F|_]=Loop, Vmap, EL, We) ->
    ?dbg("MakeEdgeLoop~n",[]), [io:format(" ~w~n", [L]) || L <- Loop],
    make_edge_loop_1(Loop, F, Vmap, EL, We);
make_edge_loop(Loop, Vmap, EL, We) ->
    %% Start with split_edge
    case lists:splitwith(fun(#{op:=Op}) -> Op =:= split_face end, Loop) of
        {FSs, []} -> split_face(FSs, Vmap, EL, We);
        {FSs, Edges} -> make_edge_loop(Edges++FSs, Vmap, EL, We)
    end.

make_edge_loop_1([#{op:=split_edge}=V1],#{op:=split_edge}=V2, Vmap, EL, We0) ->
    {We, New} = connect_verts(V1,V2,Vmap, We0),
    {[New|EL], We};
make_edge_loop_1([#{op:=split_edge}=V1|[#{op:=split_edge}=V2|_]=Rest], Last, Vmap, EL, We0) ->
    {We, New} = connect_verts(V1,V2,Vmap, We0),
    make_edge_loop_1(Rest, Last, Vmap, [New|EL], We);
make_edge_loop_1([#{op:=split_edge}=V1|Splits], Last, Vmap, EL, We0) ->
    case lists:splitwith(fun(#{op:=Op}) -> Op =:= split_face end, Splits) of
        {FSs, []} ->
            #{op:=split_edge,v:=V2} = Last,
            Face = pick_face(FSs,undefined),
            {We1, Edge} = connect_verts(V1,Last,Face,Vmap,We0),
            ok = wings_we_util:validate(We1),
            {EL0,We} = make_face_vs(FSs, array:get(V2, Vmap), Edge, Vmap, We1),
            {EL++EL0, We};
        {FSs, [#{op:=split_edge,v:=VV2}=V2|_]=Rest} ->
            Face = pick_face(FSs,undefined),
            %%?dbg("~p: ~p~n ~p~n ~p~n",[Face,E1,E2,FSs]),
            ok = wings_we_util:validate(We0),
            {We1, Edge} = connect_verts(V1,V2,Face,Vmap,We0),
            ok = wings_we_util:validate(We1),
            %%?dbg("done ~p ~p ~p ~n",[Face,V1,V2]),
            {EL0,We} = make_face_vs(FSs, array:get(VV2, Vmap), Edge, Vmap, We1),
            make_edge_loop_1(Rest, Last, Vmap, EL0++EL, We)
    end.

pick_face([#{f:=F}|Ss], undefined) ->
    pick_face(Ss, F);
pick_face([#{f:=F}|Ss], F) ->
    pick_face(Ss, F);
pick_face([], F) -> F.

connect_verts(#{v:=V1, o_n:=N1},#{v:=V2,o_n:=N2}, Vmap, We) ->
    WeV1 = array:get(V1, Vmap),
    WeV2 = array:get(V2, Vmap),
    %?dbg("~w ~w ~w ~w~n",[V1,V2,WeV1,WeV2]),
    true = is_integer(WeV1), true = is_integer(WeV2), %% Assert
    [Face|_] = [Face || {Face, [_,_]} <- wings_vertex:per_face([WeV1,WeV2],We)],
    connect_verts_1(WeV1, N1, WeV2, N2, Face, We).

connect_verts(#{v:=V1, o_n:=N1},#{v:=V2,o_n:=N2}, Face, Vmap, We) ->
    WeV1 = array:get(V1, Vmap),
    WeV2 = array:get(V2, Vmap),
    true = is_integer(WeV1), true = is_integer(WeV2), %% Assert
    %?dbg("~p ~p(~p) ~p ~p(~p) ~p ~n", [Face,V1,WeV1,N1,V2,WeV2,N2]),
    connect_verts_1(WeV1, N1, WeV2, N2, Face, We).

connect_verts_1(WeV1, N1, WeV2, N2, Face, #we{vp=Vtab}=We) ->
    case wings_vertex:edge_through(WeV1,WeV2,Face,We) of
        none when N1 =:= ignore ->
            wings_vertex:force_connect(WeV1,WeV2,Face,We);
        none ->
            N = wings_face:normal(Face, We),
            Dir = e3d_vec:cross(N,e3d_vec:sub(array:get(WeV1,Vtab),array:get(WeV2,Vtab))),
            case 0 >= e3d_vec:dot(e3d_vec:average(N1,N2), Dir) of
                true  -> wings_vertex:force_connect(WeV1,WeV2,Face,We);
                false -> wings_vertex:force_connect(WeV2,WeV1,Face,We)
            end;
        Edge ->
            %% Bugbug face should not be included in result
            {We, Edge}
    end.

make_face_vs([_]=Ss, _Vs, Edge, Vmap, We) ->
    make_face_vs_1(Ss, Edge, Vmap, [Edge], We);
make_face_vs(Ss, Vs, Edge, Vmap, #we{es=Etab}=We) ->
    case array:get(Edge, Etab) of
        #edge{vs=Vs} -> make_face_vs_1(lists:reverse(Ss), Edge, Vmap, [Edge], We);
        #edge{ve=Vs} -> make_face_vs_1(Ss, Edge, Vmap, [Edge], We)
    end.

make_face_vs_1([#{op:=split_face,v:=V}|Ss], Edge, Vmap, EL, We0) ->
    Pos = vmap_pos(V, Vmap),
    {We, New} = wings_edge:fast_cut(Edge, Pos, We0),
    make_face_vs_1(Ss, New, Vmap, [New|EL], We);
make_face_vs_1([], _, _, EL, We) ->
    {EL, We}.

split_face(Fs, Vmap, EL, We0) ->
    Face = pick_face(Fs, undefined),
    NumberOfNew = length(Fs),
    true = NumberOfNew > 2, %% Otherwise something is wrong
    We1 = wings_extrude_face:faces([Face], We0),
    FVs = wings_face:vertices_ccw(Face, We1),
    FPos = wings_face:vertex_positions(Face, We1),
    Zipped = lists:zip(FVs, FPos),
    ?dbg("Orig: ~p~n", [FVs]),
    NumberOfOld = length(FVs),
    case NumberOfOld >= NumberOfNew of
	true ->
	    KD3 = e3d_kd3:from_list(Zipped),
	    {Vs,_} = lists:mapfoldl(fun(#{v:=Vi}, Tree0) ->
					    Pos = vmap_pos(Vi, Vmap),
					    {{V,_}, Tree} = e3d_kd3:take_nearest(Pos, Tree0),
					    {{V,Pos},Tree}
				    end, KD3, Fs),
	    Vtab = lists:foldl(fun({V,Pos}, Vtab) -> array:set(V, Pos, Vtab) end,
			       We1#we.vp, Vs),
            cleanup_edges(FVs, [V||{V,_}<-Vs], Face, EL, We1#we{vp=Vtab});
	false ->
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
	    make_edge_loop(Fs1, Vmap1, EL, We1#we{vp=Vtab})
    end.

cleanup_edges(FVs, Used, Face, EL, We) ->
    %% Start with a used vertex
    {Vs1,Vs0} = lists:splitwith(fun(V) -> not lists:member(V, Used) end, FVs),
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
                   false -> {EL0,We0}
               end,
    % ?dbg("drop vs ~p~n",[Drop]),
    Es = wings_edge:from_vs(Drop, We2),
    ?dbg("drop es ~p~n",[Es]),
    We3 = wings_edge:dissolve_edges(Es, We2),
    ok = wings_we_util:validate(We3),
    {EL,We3}.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filter_tri_edges({L1,L2}) ->
    Loop = lists:zip(L1,L2),
    %?dbg("~p~n",[?FUNCTION_NAME]), [io:format("1: ~w~n2: ~w~n~n",[K1,K2])||{K1,K2}<-Loop],
    Res = filter_tri_edges_1(Loop),
    %?dbg("after ~p~n",[?FUNCTION_NAME]), [io:format("1: ~w~n2: ~w~n~n",[K1,K2])||{K1,K2}<-Res],
    lists:unzip(Res).

filter_tri_edges_1([{#{v:=V}=V1,U1}, {#{v:=V}=V2, U2}|Vs]) ->
    %% Remove edges to it self (loops)
    filter_tri_edges_1([{filter_edge(V1,V2),filter_edge(U1,U2)}|Vs]);
filter_tri_edges_1([{#{op:=split_edge,e:=none},_}|Vs]) ->
    filter_tri_edges_1(Vs);
filter_tri_edges_1([{_, #{op:=split_edge,e:=none}}|Vs]) ->
    filter_tri_edges_1(Vs);
filter_tri_edges_1([V|Vs]) ->
    [V|filter_tri_edges_1(Vs)];
filter_tri_edges_1([]) -> [].

filter_edge(_, #{op:=split_edge, e:=Edge}=V2) when Edge =/= none -> V2;
filter_edge(V1,_) -> V1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% We need to build the cycle our selfves since the edges may not be directed
%% in the correct direction. Also the V is new Vs on edges and there maybe
%% several new V on the same wings edge.

build_vtx_loops(Edges, _Acc) ->
    G = make_lookup_table(Edges),
    Comps = digraph_utils:components(G),
    ?dbg("Cs: ~w~n",[Comps]),
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
            N = if O2=:=Id -> wings_face:normal(F3, OWe);
                   true -> wings_face:normal(F2, OWe)
                end,
            Edge = wings_vertex:edge_through(A1,B1,F1,We),
            #{op=>split_edge, o=>Id, f=>F1, e=>Edge, v=>V0, vs=>Edge1, o_n=>N};
       O2 =:= Id ->
            N = if O1=:=Id -> wings_face:normal(F3, OWe);
                   true -> wings_face:normal(F1, OWe)
                end,
            Edge = wings_vertex:edge_through(A2,B2,F2,We),
            #{op=>split_edge, o=>Id, f=>F2, e=>Edge, v=>V0, vs=>Edge2, o_n=>N};
       O3 =:= Id ->
            N = wings_face:normal(F1, OWe),
            check_if_edge(#{op=>split_face, o=>Id, f=>F3, v=>V0}, N, Vmap, We)
    end;
vertex_info(#{mf1:={O1,F1}, mf2:={O2,F2}, other:={O3,F3}, p1:={{_, {A,B}=Edge0},V0}}, V0,
            Vmap, {#we{id=Id}=We,OWe}) ->
    if O1 =:= Id ->
            N = if O2=:=Id -> wings_face:normal(F3, OWe);
                   true -> wings_face:normal(F2, OWe)
                end,
            Edge = wings_vertex:edge_through(A,B,F1,We),
            #{op=>split_edge, o=>Id, f=>F1, e=>Edge, v=>V0, vs=>Edge0, o_n=>N};
       O2 =:= Id ->
            N = wings_face:normal(F1, OWe),
            check_if_edge(#{op=>split_face, o=>Id, f=>F2, v=>V0}, N, Vmap, We);
       O3 =:= Id ->
            N = wings_face:normal(F1, OWe),
            check_if_edge(#{op=>split_face, o=>Id, f=>F3, v=>V0}, N, Vmap, We)
    end;
vertex_info(#{mf2:={O1,F1}, mf1:={O2,F2}, other:={O3,F3}, p2:={{_, {A,B}=Edge0},V0}}, V0,
            Vmap, {#we{id=Id}=We,OWe}) ->
    if O1 =:= Id ->
            N = if O2=:=Id -> wings_face:normal(F3, OWe);
                   true -> wings_face:normal(F2, OWe)
                end,
            Edge = wings_vertex:edge_through(A,B,F1,We),
            #{op=>split_edge, o=>Id, f=>F1, e=>Edge, v=>V0, vs=>Edge0, o_n=>N};
       O2 =:= Id ->
            N = wings_face:normal(F1, OWe),
            check_if_edge(#{op=>split_face, o=>Id, f=>F2, v=>V0}, N, Vmap, We);
       O3 =:= Id ->
            N = wings_face:normal(F1, OWe),
            check_if_edge(#{op=>split_face, o=>Id, f=>F3, v=>V0}, N, Vmap, We)
    end.

check_if_edge(#{f:=F, v:=V}=SF, N, Vmap, #we{id=Id, vp=Vtab}=We) ->
    {Where, Pos} = array:get(V, Vmap),
    Find = fun(_,Edge,#edge{vs=V1,ve=V2},Acc) ->
                   V1P = array:get(V1, Vtab),
                   V2P = array:get(V2, Vtab),
                   case e3d_vec:line_dist_sqr(Pos, V1P, V2P) < ?EPSILON of
                       true -> [{{V1,V2}, Edge}|Acc];
                       false -> Acc
                   end
           end,
    case proplists:get_value(Id, Where) of
        undefined ->
            case wings_face:fold(Find, [], F, We) of
                [] -> SF;
                [{Vs,Edge}] -> SF#{op:=split_edge, o_n=>N, e=>Edge, vs=>Vs}
            end;
        WeV ->
            SF#{op:=split_edge, o_n=>N, e=>on_vertex, vs=>WeV}
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
    array:foldl(Add, {I0,Tree}, Vtab).

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
    Index = array:size(Vtab),
    {TessVtab, Is} = renumber_and_add_vs(Ps0, Vs, Index, Vtab, []),
    Fs = lists:foldl(fun({A,B,C}, Acc) ->
			     F = {element(A,Is), element(B, Is), element(C,Is)},
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

