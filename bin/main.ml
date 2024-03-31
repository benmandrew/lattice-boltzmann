open Lattice_boltzmann

let render_scale = 8
let height = 80
let width = 200

module D = struct
  let width = width
  let height = height
end

module BMat = Mat.BMake (D)
module Mat = Mat.FMake (D)

let viscosity = 0.02

(** Relaxation parameter *)
let omega = 1. /. ((3. *. viscosity) +. 0.5)

(** Initial and in-flow speed *)
let u0 = 0.1

let four9ths = 4. /. 9.
let one9th = 1. /. 9.
let one36th = 1. /. 36.
let n0 = Mat.init (fun _ _ -> four9ths *. (1. -. (1.5 *. (u0 ** 2.))))
let nN = Mat.init (fun _ _ -> one9th *. (1. -. (1.5 *. (u0 ** 2.))))
let nS = Mat.init (fun _ _ -> one9th *. (1. -. (1.5 *. (u0 ** 2.))))

let nE =
  Mat.init (fun _ _ -> one9th *. (1. +. (3. *. u0) +. (3. *. (u0 ** 2.))))

let nW =
  Mat.init (fun _ _ -> one9th *. (1. -. (3. *. u0) +. (3. *. (u0 ** 2.))))

let nNE =
  Mat.init (fun _ _ -> one9th *. (1. +. (3. *. u0) +. (3. *. (u0 ** 2.))))

let nSE =
  Mat.init (fun _ _ -> one9th *. (1. +. (3. *. u0) +. (3. *. (u0 ** 2.))))

let nNW =
  Mat.init (fun _ _ -> one9th *. (1. -. (3. *. u0) +. (3. *. (u0 ** 2.))))

let nSW =
  Mat.init (fun _ _ -> one9th *. (1. -. (3. *. u0) +. (3. *. (u0 ** 2.))))

let apply_barrier l lb r rb =
  Mat.init (fun i j ->
      if BMat.get lb i j && BMat.get rb i j then Mat.get r i j
      else Mat.get l i j)

(** Barrier *)
let b =
  BMat.init (fun i j ->
      i <= (height / 2) + 4
      && i >= (height / 2) - 4
      && j <= (height / 2) + 12
      && j >= (height / 2) - 12)

let bN = BMat.roll 1 H b
let bS = BMat.roll (-1) H b
let bE = BMat.roll 1 V b
let bW = BMat.roll (-1) V b
let bNE = BMat.roll 1 H b
let bNW = BMat.roll (-1) H b
let bSE = BMat.roll 1 V b
let bSW = BMat.roll (-1) V b

type state = {
  rho : Mat.t;
  ux : Mat.t;
  uy : Mat.t;
  n0 : Mat.t;
  nN : Mat.t;
  nS : Mat.t;
  nE : Mat.t;
  nW : Mat.t;
  nNE : Mat.t;
  nSE : Mat.t;
  nNW : Mat.t;
  nSW : Mat.t;
}

(** Macroscopic density *)
let get_rho { n0; nN; nS; nE; nW; nNE; nSE; nNW; nSW; _ } =
  Mat.adds [ n0; nN; nS; nE; nW; nNE; nSE; nNW; nSW ]

let get_ux { rho; nE; nW; nNE; nSE; nNW; nSW; _ } =
  let a = Mat.adds [ nE; nNE; nSE ] in
  let b = Mat.adds [ nW; nNW; nSW ] in
  Mat.(div (sub a b) rho)

let get_uy { rho; nN; nS; nNE; nSE; nNW; nSW; _ } =
  let a = Mat.adds [ nN; nNE; nNW ] in
  let b = Mat.adds [ nS; nSE; nSW ] in
  Mat.(div (sub a b) rho)

let update_macros s =
  let s = { s with rho = get_rho s } in
  let s = { s with ux = get_ux s } in
  { s with uy = get_uy s }

let state =
  update_macros
    {
      rho = Mat.zero;
      ux = Mat.zero;
      uy = Mat.zero;
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
  let nN = Mat.(roll 1 H nN) in
  let nNE = Mat.(roll 1 H nNE) in
  let nNW = Mat.(roll 1 H nNW) in
  let nS = Mat.(roll (-1) H nS) in
  let nSE = Mat.(roll (-1) H nSE) in
  let nSW = Mat.(roll (-1) H nSW) in
  let nE = Mat.(roll 1 V nE) in
  let nNE = Mat.(roll 1 V nNE) in
  let nSE = Mat.(roll 1 V nSE) in
  let nW = Mat.(roll (-1) V nW) in
  let nNW = Mat.(roll (-1) V nNW) in
  let nSW = Mat.(roll (-1) V nSW) in
  let nN = apply_barrier nN bN nS b in
  let nS = apply_barrier nS bS nN b in
  let nE = apply_barrier nE bE nW b in
  let nW = apply_barrier nW bW nE b in
  let nNE = apply_barrier nNE bNE nSW b in
  let nNW = apply_barrier nNW bNW nSE b in
  let nSE = apply_barrier nSE bSE nNW b in
  let nSW = apply_barrier nSW bSW nNE b in
  { rho; ux; uy; n0; nN; nS; nE; nW; nNE; nSE; nNW; nSW }

let collide_aux n r rho omu =
  Mat.init (fun i j ->
      ((1. -. omega) *. Mat.get n i j)
      +. (omega *. r *. Mat.get rho i j *. omu i j))

let collide s =
  let { rho; ux; uy; n0; nN; nS; nE; nW; nNE; nSE; nNW; nSW } =
    update_macros s
  in
  let ux2 = Mat.mul ux ux in
  let uy2 = Mat.mul uy uy in
  let u2 = Mat.add ux2 uy2 in
  let omu215 = Mat.init (fun i j -> 1. -. (1.5 *. Mat.get u2 i j)) in
  let uxuy = Mat.mul ux uy in
  let n0 = collide_aux n0 four9ths rho (Mat.get omu215) in
  let nN =
    collide_aux nN one9th rho (fun i j ->
        Mat.get omu215 i j +. (3. *. Mat.get uy i j) +. (4.5 *. Mat.get uy2 i j))
  in
  let nS =
    collide_aux nS one9th rho (fun i j ->
        Mat.get omu215 i j -. (3. *. Mat.get uy i j) +. (4.5 *. Mat.get uy2 i j))
  in
  let nE =
    collide_aux nE one9th rho (fun i j ->
        Mat.get omu215 i j +. (3. *. Mat.get ux i j) +. (4.5 *. Mat.get ux2 i j))
  in
  let nW =
    collide_aux nW one9th rho (fun i j ->
        Mat.get omu215 i j -. (3. *. Mat.get ux i j) +. (4.5 *. Mat.get ux2 i j))
  in
  let nNE =
    collide_aux nNE one36th rho (fun i j ->
        Mat.get omu215 i j
        +. (3. *. (Mat.get ux i j +. Mat.get uy i j))
        +. (4.5 *. Mat.get u2 i j)
        +. (9. *. Mat.get uxuy i j))
  in
  let nNW =
    collide_aux nNW one36th rho (fun i j ->
        Mat.get omu215 i j
        +. (3. *. (Mat.get uy i j -. Mat.get ux i j))
        +. (4.5 *. Mat.get u2 i j)
        -. (9. *. Mat.get uxuy i j))
  in
  let nSE =
    collide_aux nSE one36th rho (fun i j ->
        Mat.get omu215 i j
        +. (3. *. (Mat.get ux i j -. Mat.get uy i j))
        +. (4.5 *. Mat.get u2 i j)
        -. (9. *. Mat.get uxuy i j))
  in
  let nSW =
    collide_aux nSW one36th rho (fun i j ->
        Mat.get omu215 i j
        -. (3. *. (Mat.get ux i j +. Mat.get uy i j))
        +. (4.5 *. Mat.get u2 i j)
        +. (9. *. Mat.get uxuy i j))
  in
  (* Steady rightward flow *)
  let nE =
    Mat.init (fun i j ->
        if i = 0 then
          1. +. (3. *. u0) +. (4.5 *. (u0 ** 2.)) -. (1.5 *. (u0 ** 2.))
        else Mat.get nE i j)
  in
  let nW =
    Mat.init (fun i j ->
        if i = 0 then
          one9th
          *. (1. -. (3. *. u0) +. (4.5 *. (u0 ** 2.)) -. (1.5 *. (u0 ** 2.)))
        else Mat.get nW i j)
  in
  let nNE =
    Mat.init (fun i j ->
        if i = 0 then
          one36th
          *. (1. +. (3. *. u0) +. (4.5 *. (u0 ** 2.)) -. (1.5 *. (u0 ** 2.)))
        else Mat.get nNE i j)
  in
  let nSE =
    Mat.init (fun i j ->
        if i = 0 then
          one36th
          *. (1. +. (3. *. u0) +. (4.5 *. (u0 ** 2.)) -. (1.5 *. (u0 ** 2.)))
        else Mat.get nSE i j)
  in
  let nNW =
    Mat.init (fun i j ->
        if i = 0 then
          one36th
          *. (1. -. (3. *. u0) +. (4.5 *. (u0 ** 2.)) -. (1.5 *. (u0 ** 2.)))
        else Mat.get nNW i j)
  in
  let nSW =
    Mat.init (fun i j ->
        if i = 0 then
          one36th
          *. (1. -. (3. *. u0) +. (4.5 *. (u0 ** 2.)) -. (1.5 *. (u0 ** 2.)))
        else Mat.get nSW i j)
  in
  { rho; ux; uy; n0; nN; nS; nE; nW; nNE; nSE; nNW; nSW }

let curl ux uy =
  let a = Mat.(add (roll (-1) V uy) (roll 1 V uy)) in
  let b = Mat.(add (roll (-1) H ux) (roll 1 H ux)) in
  Mat.sub a b

let next_frame s = collide @@ stream s

module G = Graphics

let exit_handler status =
  if status.G.keypressed && status.key == ' ' then raise Exit else ()

let discretise scale m =
  let lum =
    Array.init height (fun j ->
        Array.init width (fun i ->
            int_of_float @@ Float.mul 255. @@ Mat.get m i j))
  in
  Array.init (height * scale) (fun j ->
      let j' = j / scale in
      Array.init (width * scale) (fun i ->
          let i' = i / scale in
          let barrier = if BMat.get b i' j' then 120 else 0 in
          let lum = Array.get (Array.get lum j') i' in
          if lum >= 0 then G.rgb lum 0 barrier else G.rgb 0 (abs lum) barrier))

let rec loop s =
  G.draw_image (G.make_image @@ discretise render_scale @@ curl s.ux s.uy) 0 0;
  loop @@ next_frame s

let () =
  ignore (state, stream, collide, curl, discretise, exit_handler);
  G.open_graph "";
  G.resize_window (width * render_scale) (height * render_scale);
  loop state
(* G.loop_at_exit [ Key_pressed ] exit_handler *)
