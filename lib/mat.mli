module type D = sig
  val width : int
  val height : int
end

module FMake : functor (_ : D) -> sig
  type t

  val init : (int -> int -> float) -> t
  val zero : t
  val get : t -> int -> int -> float
  val op : (float -> float -> float) -> t -> t -> t
  val add : t -> t -> t
  val sub : t -> t -> t
  val mul : t -> t -> t
  val div : t -> t -> t
  val adds : t list -> t

  type dir = H | V

  val modulo : int -> int -> int
  val roll : int -> dir -> t -> t
end

module BMake : functor (_ : D) -> sig
  type t

  val init : (int -> int -> bool) -> t
  val zero : t
  val get : t -> int -> int -> bool
  val op : (bool -> bool -> bool) -> t -> t -> t

  type dir = H | V

  val modulo : int -> int -> int
  val roll : int -> dir -> t -> t
end
