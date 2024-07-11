load("test_pairing.sage")
load("global_parameters.sage")
load("r1cs_qap.sage")


def setup(k=2, apoly, bpoly, cpoly, n):
    random_variables = [int(FF.random_element()) for _ in range(5)]
    tau, alpha, beta, gamma, delta = random_variables

    k_vk, k_pk, z_d = [], [], []
    apoly1, bpoly1, bpoly2 = [], [], []
    for i in range(0, n+1):
        apoly1.append(multiply(G1, Pr_k(apoly[i].list())(tau)))
        bpoly1.append(multiply(G1, Pr_k(bpoly[i].list())(tau)))
        bpoly2.append(multiply(G2, Pr_k(bpoly[i].list())(tau)))
        a_it = beta * Pr_k(apoly[i].list())(tau)
        b_it = alpha * Pr_k(bpoly[i].list())(tau)
        c_it = Pr_k(cpoly[i].list())(tau)
        if i <= k:
            k_vki = 1/gamma * (beta * a_it + alpha * b_it + c_it)
            k_vk.append(multiply(G1, k_vki))
        else:
            k_pki = 1/delta * (beta * a_it + alpha * b_it + c_it)
            k_pk.append(multiply(G1, k_pki))

    for j in range(0, m-2):
        z_d.append(multiply(G1, 1/delta * t^j * Z_d(tau)))
    prover_key = []
    verifier_key = []

    alpha1 = multply(G1, alpha)
    beta1 = multiply(G1, beta)
    beta2 = multiply(G2, b2)
    delta1 = multiply(G1, d)
    delta2 = multiply(G2, d)
    gamma2 = multiply(G2, gamma)

    prover_key.append([alpha1, beta1, beta2, delta1, delta2])
    prover_key.append(apoly1)
    prover_key.append(bpoly1)
    prover_key.append(bpoly2)
    prover_key.append(k_pk)
    prover_key.append(z_d)

    verifier_key.append(pairing(a1, b2))
    verifier_key.append(gamma2)
    verifier_key.append(delta2)
    verifier_key.append(k_vk)

    return prover_key, verifier_key


def prover(prover_key, z, h):
    r, s = FF.random_element(), FF.random_element()
    alpha1, beta1, beta2, delta1, delta2 = prover_key[0]
    apoly1 = prover_key[1]
    bpoly1 = prover_key[2]
    bpoly2 = prover_key[3]
    k_pk = prover_key[4]
    z_d = prover_key[5]

    rd1, rd2 = multiply(delta1, r)
    sd1, sd2 = multiply(delta1, s), multiply(delta2, s)
    ar1 = add(alpha1, rd1)
    bs1 = add(beta1, sd1)
    bs2 = add(beta2, sd2)
    for i in range(len(apoly1)):
        ar1 = add(ar1, multiply(z[i], apoly1[i]))
        bs1 = add(bs1, multiply(z[i], bpoly1[i]))
        bs2 = add(bs2, multiply(z[i], bpoly2[i]))

    krs1 = multiply(ar1, s)
    krs1 = add(krs1, multiply(bs1, r))
    rsd1 = multiply(sd1, r)
    krs1 = add(krs1, neg(rsd1))

    for i in range(len(k_pk)):
        krs1 = add(krs1, multiply(k_pk[i], z[i]))
    for i in range(len(z_d)):
        krs1 = add(krs1, multiply(z_d[i], h[i]))

    return [ar1, bs2, krs1]


def verifier(verifier_key, proof, x):
    ab_pairing, gamma2, delta2, k_vk = verifier_key
    ar1, bs2, krs1 = proof
    p1 = pairing(ar1, bs2)
    p2 = ab_pairing

    kx_sum = multiply(k_vk[0], x[0])
    for i in range(1, len(k_vk)):
        kx_sum = add(kx_sum, multiply(k_vk[i], x[i]))
    p2 += pairing(kx_sum, gamma2)
    p2 += pairing(krs1, delta2)
    return p1 == p2


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

    Lp = matrix_to_polynomials(L, FF)
    Rp = matrix_to_polynomials(R, FF)
    Op = matrix_to_polynomials(O, FF)

    A = PR_k((w * Lp).list())
    B = PR_k((w * Rp).list())
    C = PR_k((w * Op).list())
    target_poly, h = create_vanishing_polynomial(A, B, C, L.nrows())



if __name__ == "__main__":
    main()
