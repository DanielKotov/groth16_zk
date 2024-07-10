def matrix_to_polynomials(matrix, FF):
    row, col = matrix.nrows(), matrix.ncols()
    poly_list = []
    for i in range(col):
        points = []
        for j in range(row):
            x = FF(j+1)
            y = matrix[j][i]
            points.append((x, y))
        poly = PR_k.lagrange_polynomial(points)
        coefs = poly.list()[::-1]
        print(f"coefs == {coefs}")
        if len(coefs) < row:
            coefs = [0] * (row - len(coefs)) + coefs
        print(coefs)
        poly_list.append(coefs[::-1])
    return Matrix(FF, poly_list)



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
        tau = FF.random_element()
        assert(U * V == W + h * T)
    print("GOOOOL")
