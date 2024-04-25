open Lattice_boltzmann
(* open Lattice_boltzmann *)

let render_scale = 16
let height = 40
let width = 100

module D = struct
  let width = width
  let height = height
end

module Sim = Simulate.Make (D)
module FMat = Sim.FMat
module BMat = Sim.BMat

let next_frame s = Sim.(collide @@ update_macroscopics @@ stream s)

module G = Graphics
module A = Array

let exit_handler status =
  if status.G.keypressed && status.key == ' ' then raise Exit else ()

let discretise scale m =
  let f x =
    let x = 8. *. x in
    if x < 0. then -.x /. (x -. 1.) else x /. (x +. 1.)
  in
  let lum =
    A.init height (fun j ->
        A.init width (fun i ->
            int_of_float @@ Float.mul 255. @@ f @@ FMat.get m i j))
  in
  A.init (height * scale) (fun j ->
      let j' = j / scale in
      A.init (width * scale) (fun i ->
          let i' = i / scale in
          let barrier = if BMat.get Sim.b i' j' then 120 else 0 in
          let lum = A.get (A.get lum j') i' in
          if lum >= 0 then G.rgb lum 0 barrier else G.rgb 0 (abs lum) barrier))

let print_variance i s c =
  let sum = FMat.fold ( +. ) 0. s.Sim.rho in
  let n = float_of_int (width * height) in
  let mean = sum /. n in
  let variance =
    FMat.fold (fun acc x -> acc +. ((x -. mean) ** 2.)) 0. c /. n
  in
  Printf.printf "\rframe %d - VAR(curl) = %f" i variance

let rec loop i s =
  let c = Sim.curl s.Sim.ux s.uy in
  G.draw_image (G.make_image @@ discretise render_scale c) 0 0;
  if i mod 10 == 0 then (
    print_variance i s c;
    flush stdout);
  loop (i + 1) @@ next_frame s

let () =
  ignore (exit_handler);
  G.open_graph "";
  G.resize_window (width * render_scale) (height * render_scale);
  loop 0 Sim.state
(* G.loop_at_exit [ Key_pressed ] exit_handler *)
