module type V = sig
  type t
end

module type D = sig
  val width : int
  val height : int
end

module Make (V : V) (D : D) = struct
  type t = V.t array array

  let init v =
    Array.init D.height (fun j -> Array.init D.width (fun i -> v i j))

  let get m i j = Array.(get (get m j) i)
  let op f = Array.map2 (fun a b -> Array.map2 (fun a b -> f a b) a b)

  type dir = H | V

  let modulo x y =
    let r = x mod y in
    if r >= 0 then r else r + y

  let roll n d m =
    let f i j =
      match d with
      | H -> get m (modulo (i + n) D.width) j
      | V -> get m i (modulo (j + n) D.height)
    in
    init f
end

module FMake (D : D) = struct
  include Make (Float) (D)

  let zero = init (fun _ _ -> 0.)
  let add = op Float.add
  let sub = op Float.sub
  let mul = op Float.mul
  let div = op Float.div
  let adds = List.fold_left (fun acc m -> add acc m) zero
end

module BMake (D : D) = struct
  include Make (Bool) (D)

  let zero = init (fun _ _ -> false)
end
