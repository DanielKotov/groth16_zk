def matrix_to_polynomials(matrix):
    row, col = matrix.nrows(), matrix.ncols()
    poly_list = []
    for i in range(col):
        points = []
        for j in range(row):
            x = F(j+1)
            y = matrix[j][i]
            points.append((x, y))
        poly = P.lagrange_polynomial(points)
        coefs = poly.list()[::-1]
        print(f"coefs == {coefs}")
        if len(coefs) < row:
            coefs = [0] * (row - len(coefs)) + coefs
        print(coefs)
        poly_list.append(coefs[::-1])
    return Matrix(F, poly_list)



def create_vanishing_polynomial(U, V, W, matrix):
    T = (k - 1)
    row = matrix.nrows()
    for j in range(2, row+1):
        T *= (k - j)
    print("\nT=", T)
    for i in range(1, row+1):
        print(f"T({i}) == ", T(i))
    h = (U * V - W) // T
    rem = (U * V - W) % T
    return T, h, rem


def test_vanishing_polynomial(U, V, W, matrix):
    T, h, rem = create_vanishing_polynomial(U, V, W, matrix)
    for i in range(1, 100):
        tau = F.random_element()
        assert(U * V == W + h * T)
    print("GOOOOL")


if __name__ == "__main__":

    p = 71
    F = GF(p)
    P.<k> = PolynomialRing(F)
    L = Matrix(F, [[0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 5, 0, 0, 0, 0, 0], [0, 0, 0, 0, 4, 0, 0, 0], [0, 0, 13, 0, 0, 0, 0, 0]])
    R = Matrix(F, [[0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0]])
    O = Matrix(F, [[0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1], [0, 1, 0, 10, p - 1, 0, p - 1, 1]])

    x = F(2)
    y = F(3)
    v1 = x * x
    v2 = y * y
    v3 = 5 * x * v1
    v4 = 4 * v1 * v2
    out = 5*x**3 - 4*x**2*y**2 + 13*x*y**2 + x**2 - 10*y
    w = vector(F, [1, out, x, y, v1, v2, v3, v4])

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
    U = P((w * Lp).list())
    V = P((w * Rp).list())
    W = P((w * Op).list())

    print(U)
    print(V)
    print(W)

    test_vanishing_polynomial(U, V, W, L)
