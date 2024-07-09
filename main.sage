load("test_pairing.sage")
load("global_parameters.sage")
load("r1cs_qap.sage")

def main():
    L = Matrix(FF, [[0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 5, 0, 0, 0, 0, 0], [0, 0, 0, 0, 4, 0, 0, 0], [0, 0, 13, 0, 0, 0, 0, 0]])
    R = Matrix(FF, [[0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0]])
    O = Matrix(FF, [[0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1], [0, 1, 0, 10, p - 1, 0, p - 1, 1]])

    x = FF(2)
    y = FF(3)
    v1 = x * x
    v2 = y * y
    v3 = 5 * x * v1
    v4 = 4 * v1 * v2
    out = 5*x**3 - 4*x**2*y**2 + 13*x*y**2 + x**2 - 10*y
    w = vector(FF, [1, out, x, y, v1, v2, v3, v4])

    Lw = L * w
    Rw = R * w
    Ow = O * w

    LRw = Lw.pairwise_product(Rw)

    print(f"w = {w}")

    print(f"Lw = {Lw}")
    print(f"Rw = {Rw}")
    print(f"Ow = {Ow}")

    print(f"LRw = {LRw}")
    print(f"Ow = {Ow}")

    Lp = matrix_to_polynomials(L)
    Rp = matrix_to_polynomials(R)
    Op = matrix_to_polynomials(O)
    print()
    print(Lp)
    print()
    print(Rp)
    print()
    print(Op)
    print()
    U = PR((w * Lp).list())
    V = PR((w * Rp).list())
    W = PR((w * Op).list())
    print(U)
    print(V)
    print(W)
    test_vanishing_polynomial(U, V, W, L)


if __name__ == "__main__":
    main()
