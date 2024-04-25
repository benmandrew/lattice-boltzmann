module Make (D : S.D) = struct
  module FMat = Mat.FMake (D)
  module BMat = Mat.BMake (D)

  let viscosity = 0.02

  (** Relaxation parameter *)
  let omega = 1. /. ((3. *. viscosity) +. 0.5)

  (** Initial and in-flow speed *)
  let u0 = -0.1

  let four9ths = 4. /. 9.
  let one9th = 1. /. 9.
  let one36th = 1. /. 36.
  let n0 = FMat.init (fun _ _ -> four9ths *. (1. -. (1.5 *. (u0 ** 2.))))
  let nN = FMat.init (fun _ _ -> one9th *. (1. -. (1.5 *. (u0 ** 2.))))
  let nS = FMat.init (fun _ _ -> one9th *. (1. -. (1.5 *. (u0 ** 2.))))

  let nE =
    FMat.init (fun _ _ -> one9th *. (1. +. (3. *. u0) +. (3. *. (u0 ** 2.))))

  let nW =
    FMat.init (fun _ _ -> one9th *. (1. -. (3. *. u0) +. (3. *. (u0 ** 2.))))

  let nNE =
    FMat.init (fun _ _ -> one36th *. (1. +. (3. *. u0) +. (3. *. (u0 ** 2.))))

  let nSE =
    FMat.init (fun _ _ -> one36th *. (1. +. (3. *. u0) +. (3. *. (u0 ** 2.))))

  let nNW =
    FMat.init (fun _ _ -> one36th *. (1. -. (3. *. u0) +. (3. *. (u0 ** 2.))))

  let nSW =
    FMat.init (fun _ _ -> one36th *. (1. -. (3. *. u0) +. (3. *. (u0 ** 2.))))

  let apply_barrier l lb r rb =
    FMat.init (fun i j ->
        if BMat.get lb i j && BMat.get rb i j then FMat.get r i j
        else FMat.get l i j)

  (** Barrier *)
  let b =
    BMat.init (fun i j ->
        i <= (D.height / 2) + 4
        && i >= (D.height / 2) - 4
        && j <= (D.height / 2) + 12
        && j >= (D.height / 2) - 12)

  let bN = BMat.roll 1 V b
  let bS = BMat.roll (-1) V b
  let bE = BMat.roll 1 H b
  let bW = BMat.roll (-1) H b
  let bNE = BMat.roll 1 H bN
  let bNW = BMat.roll (-1) H bN
  let bSE = BMat.roll 1 H bS
  let bSW = BMat.roll (-1) H bS

  type state = {
    rho : FMat.t;
    ux : FMat.t;
    uy : FMat.t;
    n0 : FMat.t;
    nN : FMat.t;
    nS : FMat.t;
    nE : FMat.t;
    nW : FMat.t;
    nNE : FMat.t;
    nSE : FMat.t;
    nNW : FMat.t;
    nSW : FMat.t;
  }

  (** Macroscopic density *)
  let get_rho { n0; nN; nS; nE; nW; nNE; nSE; nNW; nSW; _ } =
    FMat.adds [ n0; nN; nS; nE; nW; nNE; nSE; nNW; nSW ]

  let get_ux { rho; nE; nW; nNE; nSE; nNW; nSW; _ } =
    let a = FMat.adds [ nE; nNE; nSE ] in
    let b = FMat.adds [ nW; nNW; nSW ] in
    FMat.(div (sub a b) rho)

  let get_uy { rho; nN; nS; nNE; nSE; nNW; nSW; _ } =
    let a = FMat.adds [ nN; nNE; nNW ] in
    let b = FMat.adds [ nS; nSE; nSW ] in
    FMat.(div (sub a b) rho)

  let update_macroscopics s =
    let s = { s with rho = get_rho s } in
    let s = { s with ux = get_ux s } in
    { s with uy = get_uy s }

  let state =
    {
      rho = FMat.zero;
      ux = FMat.zero;
      uy = FMat.zero;
      n0;
      nN;
      nS;
      nE;
      nW;
      nNE;
      nSE;
      nNW;
      nSW;
    }

  let stream { rho; ux; uy; n0; nN; nS; nE; nW; nNE; nSE; nNW; nSW } =
    let nN = FMat.(roll 1 V nN) in
    let nNE = FMat.(roll 1 V nNE) in
    let nNW = FMat.(roll 1 V nNW) in
    let nS = FMat.(roll (-1) V nS) in
    let nSE = FMat.(roll (-1) V nSE) in
    let nSW = FMat.(roll (-1) V nSW) in
    let nE = FMat.(roll 1 H nE) in
    let nNE = FMat.(roll 1 H nNE) in
    let nSE = FMat.(roll 1 H nSE) in
    let nW = FMat.(roll (-1) H nW) in
    let nNW = FMat.(roll (-1) H nNW) in
    let nSW = FMat.(roll (-1) H nSW) in
    let nN = apply_barrier nN bN nS b in
    let nS = apply_barrier nS bS nN b in
    let nE = apply_barrier nE bE nW b in
    let nW = apply_barrier nW bW nE b in
    let nNE = apply_barrier nNE bNE nSW b in
    let nNW = apply_barrier nNW bNW nSE b in
    let nSE = apply_barrier nSE bSE nNW b in
    let nSW = apply_barrier nSW bSW nNE b in
    { rho; ux; uy; n0; nN; nS; nE; nW; nNE; nSE; nNW; nSW }

  (** Steady rightward flow *)
  let create_flow { rho; ux; uy; n0; nN; nS; nE; nW; nNE; nSE; nNW; nSW } =
    let nE =
      FMat.init (fun i j ->
          if i = 0 then
            one9th
            *. (1. +. (3. *. u0) +. (4.5 *. (u0 ** 2.)) -. (1.5 *. (u0 ** 2.)))
          else FMat.get nE i j)
    in
    let nW =
      FMat.init (fun i j ->
          if i = 0 then
            one9th
            *. (1. -. (3. *. u0) +. (4.5 *. (u0 ** 2.)) -. (1.5 *. (u0 ** 2.)))
          else FMat.get nW i j)
    in
    let nNE =
      FMat.init (fun i j ->
          if i = 0 then
            one36th
            *. (1. +. (3. *. u0) +. (4.5 *. (u0 ** 2.)) -. (1.5 *. (u0 ** 2.)))
          else FMat.get nNE i j)
    in
    let nSE =
      FMat.init (fun i j ->
          if i = 0 then
            one36th
            *. (1. +. (3. *. u0) +. (4.5 *. (u0 ** 2.)) -. (1.5 *. (u0 ** 2.)))
          else FMat.get nSE i j)
    in
    let nNW =
      FMat.init (fun i j ->
          if i = 0 then
            one36th
            *. (1. -. (3. *. u0) +. (4.5 *. (u0 ** 2.)) -. (1.5 *. (u0 ** 2.)))
          else FMat.get nNW i j)
    in
    let nSW =
      FMat.init (fun i j ->
          if i = 0 then
            one36th
            *. (1. -. (3. *. u0) +. (4.5 *. (u0 ** 2.)) -. (1.5 *. (u0 ** 2.)))
          else FMat.get nSW i j)
    in
    { rho; ux; uy; n0; nN; nS; nE; nW; nNE; nSE; nNW; nSW }

  let collide_aux n r rho omu =
    FMat.init (fun i j ->
        ((1. -. omega) *. FMat.get n i j)
        +. (omega *. r *. FMat.get rho i j *. omu i j))

  let collide { rho; ux; uy; n0; nN; nS; nE; nW; nNE; nSE; nNW; nSW } =
    let ux2 = FMat.mul ux ux in
    let uy2 = FMat.mul uy uy in
    let u2 = FMat.add ux2 uy2 in
    let omu215 = FMat.init (fun i j -> 1. -. (1.5 *. FMat.get u2 i j)) in
    let uxuy = FMat.mul ux uy in
    let n0 = collide_aux n0 four9ths rho (FMat.get omu215) in
    let nN =
      collide_aux nN one9th rho (fun i j ->
          FMat.get omu215 i j
          +. (3. *. FMat.get uy i j)
          +. (4.5 *. FMat.get uy2 i j))
    in
    let nS =
      collide_aux nS one9th rho (fun i j ->
          FMat.get omu215 i j
          -. (3. *. FMat.get uy i j)
          +. (4.5 *. FMat.get uy2 i j))
    in
    let nE =
      collide_aux nE one9th rho (fun i j ->
          FMat.get omu215 i j
          +. (3. *. FMat.get ux i j)
          +. (4.5 *. FMat.get ux2 i j))
    in
    let nW =
      collide_aux nW one9th rho (fun i j ->
          FMat.get omu215 i j
          -. (3. *. FMat.get ux i j)
          +. (4.5 *. FMat.get ux2 i j))
    in
    let nNE =
      collide_aux nNE one36th rho (fun i j ->
          FMat.get omu215 i j
          +. (3. *. (FMat.get ux i j +. FMat.get uy i j))
          +. (4.5 *. FMat.get u2 i j)
          +. (9. *. FMat.get uxuy i j))
    in
    let nNW =
      collide_aux nNW one36th rho (fun i j ->
          FMat.get omu215 i j
          +. (3. *. (FMat.get uy i j -. FMat.get ux i j))
          +. (4.5 *. FMat.get u2 i j)
          -. (9. *. FMat.get uxuy i j))
    in
    let nSE =
      collide_aux nSE one36th rho (fun i j ->
          FMat.get omu215 i j
          +. (3. *. (FMat.get ux i j -. FMat.get uy i j))
          +. (4.5 *. FMat.get u2 i j)
          -. (9. *. FMat.get uxuy i j))
    in
    let nSW =
      collide_aux nSW one36th rho (fun i j ->
          FMat.get omu215 i j
          -. (3. *. (FMat.get ux i j +. FMat.get uy i j))
          +. (4.5 *. FMat.get u2 i j)
          +. (9. *. FMat.get uxuy i j))
    in
    create_flow { rho; ux; uy; n0; nN; nS; nE; nW; nNE; nSE; nNW; nSW }


    let curl ux uy =
      let a = FMat.(add (roll (-1) V uy) (roll 1 H ux)) in
      let b = FMat.(add (roll 1 V uy) (roll (-1) H ux)) in
      FMat.sub a b

end
