from py_ecc.optimized_bn128 import multiply, G1, G2, add, pairing, neg, normalize


def commit_polynomial_in_point(f, tau, g):
# want to commit evaluation of polynomial in the point tau using the next formula:
# cm(tau) = a_0 * G_0 + a_1 * G_1 + ... + a_n * P_n,
# where P_i = tau^i * G, so, cm(tau) = f(tau) * G
    coeff = [int(elem) for elem in f.list()]
    d = f.degree()
    points = []
    for i in range(0, d):
        point = multiply(g, int(tau^(i+1)))
        points.append(multiply(point, coeff[i]))
    res = points[0]
    for i in range(1, len(points)):
        res = add(res, points[i])
    return res
