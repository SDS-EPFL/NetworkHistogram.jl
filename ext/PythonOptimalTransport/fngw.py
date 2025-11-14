import numpy as np

from ot.utils import dist, UndefinedParameter, list_to_array
from ot.utils import check_random_state
from ot.backend import get_backend
from ot.optim import line_search_armijo, solve_1d_linesearch_quad
from ot.lp import emd

from ot.gromov import init_matrix as init_matrix_A
from ot.gromov import gwloss, gwggrad


def fngw_barycenters(
    N,
    Fs,
    As,
    Cs,
    ps,
    lambdas,
    alpha,
    beta,
    fixed_structure=False,
    fixed_node_features=False,
    fixed_edge_features=False,
    p=None,
    dist_fun_C="l2_norm",
    dist_fun_A="square_loss",
    max_iter=100,
    tol=1e-9,
    verbose=False,
    log=False,
    init_C=None,
    init_F=None,
    init_A=None,
    random_state=None,
):
    r"""Compute the FNGW barycenter as presented eq (12) in our paper

    Parameters
    ----------
    N : int
        Desired number of samples of the target barycenter
    Fs: list of array-like, each element has shape (ns,d)
        Node features of all samples
    As : list of array-like, each element has shape (ns,ns)
        Structure matrices of all samples
    Cs : list of array-like, each element has shape (ns,ns,d')
        Edge feature tensors of all samples
    ps : list of array-like, each element has shape (ns,)
        Masses of all samples.
    lambdas : list of float
        List of the `S` spaces' weights
    alpha : float
        Alpha parameter for the FNGW distance
    beta : float
        Alpha parameter for the FNGW distance
    fixed_structure : bool
        Whether to fix the structure of the barycenter during the updates
    fixed_node_features : bool
        Whether to fix the node feature of the barycenter during the updates
    fixed_edge_features : bool
        Whether to fix the edge feature of the barycenter during the updates
    dist_fun_A : str
        Loss function used for the solver either 'square_loss' or 'kl_loss'
    dist_fun_C : str
        Inner distance function used for the solver. Now only 'l2_norm'
    max_iter : int, optional
        Max number of iterations
    tol : float, optional
        Stop threshold on error (>0).
    verbose : bool, optional
        Print information along iterations.
    log : bool, optional
        Record log if True.
    init_C : array-like, shape (N,N,d'), optional
        Initialization for the barycenters' edge feature tensor. If not set
        a random init is used.
    init_F : array-like, shape (N,d), optional
        Initialization for the barycenters' node features. If not set a
        random init is used.
    init_A : array-like, shape (N,N), optional
        Initialization for the barycenters' structure matrix. If not set
        a random init is used.
    random_state : int or RandomState instance, optional
        Fix the seed for reproducibility

    Returns
    -------
    F : array-like, shape (`N`, `d`)
        Barycenters' features
    A : array-like, shape (`N`, `N`)
        Barycenters' structure matrix
    C : array-like, shape (`N`, `N`, `d'`)
        Barycenters' edge feature tensor
    log : dict
        Only returned when log=True. It contains the keys:

        - :math:`\mathbf{T}`: list of (`N`, `ns`) transport matrices
        - :math:`(\mathbf{M}_s)_s`: all distance matrices between the feature
        of the barycenter and the other features
        :math:`(dist(\mathbf{X}, \mathbf{Y}_s))_s` shape (`N`, `ns`)

    """
    Cs = list_to_array(*Cs)
    As = list_to_array(*As)
    ps = list_to_array(
        *ps
    )  # list to array bug when only one list has length one
    Fs = list_to_array(*Fs)
    if not isinstance(Cs, list):
        Cs = [Cs]
    if not isinstance(As, list):
        As = [As]
    if not isinstance(ps, list):
        ps = [ps]
    if not isinstance(Fs, list):
        Fs = [Fs]

    p = list_to_array(p)
    nx = get_backend(*Cs, *Fs, *As, *ps)

    S = len(Cs)
    d = Fs[0].shape[1]  # dimension on the node features
    d_edge = Cs[0].shape[2]
    if p is None:
        p = nx.ones(N, type_as=Cs[0]) / N

    if fixed_edge_features:
        if init_C is None:
            raise UndefinedParameter("If C is fixed it must be initialized")
        else:
            C = init_C
    else:
        if init_C is None:
            generator = check_random_state(random_state)
            C = generator.randn(N, N, d_edge)
            C = nx.from_numpy(C, type_as=ps[0])
        else:
            C = init_C

    if fixed_structure:
        if init_A is None:
            raise UndefinedParameter("If A is fixed it must be initialized")
        else:
            A = init_A
    else:
        if init_A is None:
            generator = check_random_state(random_state)
            xalea = generator.randn(N, 2)
            A = dist(xalea, xalea)
            A = nx.from_numpy(A, type_as=ps[0])
        else:
            A = init_A

    if fixed_node_features:
        if init_F is None:
            raise UndefinedParameter("If F is fixed it must be initialized")
        else:
            F = init_F
    else:
        if init_F is None:
            F = nx.zeros((N, d), type_as=ps[0])
        else:
            F = init_F

    T = [nx.outer(p, q) for q in ps]

    Ms = [dist(F, Fs[s]) for s in range(len(Fs))]

    cpt = 0
    err_node_feature = 1
    err_structure = 1
    err_edge_feature = 1

    if log:
        log_ = {}
        log_["err_node_feature"] = []
        log_["err_edge_feature"] = []
        log_["err_structure"] = []
        log_["Ts_iter"] = []

    while (
        err_node_feature > tol or err_structure > tol or err_edge_feature > tol
    ) and cpt < max_iter:
        Cprev = C
        Aprev = A
        Xprev = F

        if not fixed_node_features:
            Fs_temp = [y.T for y in Fs]
            F = update_node_feature_matrix(lambdas, Fs_temp, T, p).T

        Ms = [dist(F, Fs[s]) for s in range(len(Fs))]

        if not fixed_structure:
            if dist_fun_A == "square_loss":
                T_temp = [t.T for t in T]
                A = update_structure_matrix(p, lambdas, T_temp, As)

        if not fixed_edge_features:
            if dist_fun_C == "l2_norm":
                T_temp = [t.T for t in T]
                C = update_edge_feature_tensor(p, lambdas, T_temp, Cs)

        T = [
            fused_network_gromov_wasserstein2(
                Ms[s],
                C,
                Cs[s],
                A,
                As[s],
                p,
                ps[s],
                dist_fun_C,
                dist_fun_A,
                alpha,
                beta,
                numItermax=max_iter,
                stopThr=1e-5,
                verbose=verbose > 2,
                log=True,
            )[1]["T"]
            for s in range(S)
        ]

        # T is N,ns
        err_node_feature = nx.norm(F - nx.reshape(Xprev, (N, d)))
        err_structure = nx.norm(A - Aprev)
        err_edge_feature = nx.norm(C - Cprev)
        if log:
            log_["err_node_feature"].append(err_node_feature)
            log_["err_edge_feature"].append(err_edge_feature)
            log_["err_structure"].append(err_structure)
            log_["Ts_iter"].append(T)

        if verbose:
            if cpt % 200 == 0:
                print("{:5s}|{:12s}".format("It.", "Err") + "\n" + "-" * 19)
            print("{:5d}|{:8e}|".format(cpt, err_structure))
            print("{:5d}|{:8e}|".format(cpt, err_node_feature))
            print("{:5d}|{:8e}|".format(cpt, err_edge_feature))
            print("\n")

        cpt += 1

    if log:
        log_["T"] = T  # from target to Fs
        log_["p"] = p
        log_["Ms"] = Ms

    if log:
        return F, A, C, log_
    else:
        return F, A, C


def update_structure_matrix(p, lambdas, T, As):
    r"""Updates :math:`\mathbf{C}` according to the L2 Loss kernel with the
    `S` :math:`\mathbf{T}_s` couplings.
    It is calculated at each iteration
    Parameters
    ----------
    p : array-like, shape (N,)
        Masses in the targeted barycenter.
    lambdas : list of float
        List of the `S` spaces' weights.
    T : list of S array-like of shape (ns, N)
        The `S` :math:`\mathbf{T}_s` couplings calculated at each iteration.
    As : list of S array-like, shape (ns, ns)
        Metric cost matrices.
    Returns
    -------
    A : array-like, shape (`nt`, `nt`)
        Updated :math:`\mathbf{A}` matrix.
    """
    p = list_to_array(p)
    T = list_to_array(*T)
    As = list_to_array(*As)
    nx = get_backend(*As, *T, p)

    tmpsum = sum(
        [
            lambdas[s] * nx.dot(nx.dot(T[s].T, As[s]), T[s])
            for s in range(len(T))
        ]
    )
    ppt = nx.outer(p, p)
    return tmpsum / ppt


def update_node_feature_matrix(lambdas, Fs, Ts, p):
    r"""Updates the feature with respect to the `S` :math:`\mathbf{T}_s`
    couplings.
    Parameters
    ----------
    p : array-like, shape (N,)
        masses in the targeted barycenter
    lambdas : list of float
        List of the `S` spaces' weights
    Ts : list of S array-like, shape (ns,N)
        The `S` :math:`\mathbf{T}_s` couplings calculated at each iteration
    Fs : list of S array-like, shape (d,ns)
        The features.

    Returns
    -------
    F : array-like, shape (`d`, `N`)
    """
    p = list_to_array(p)
    Ts = list_to_array(*Ts)
    Fs = list_to_array(*Fs)
    if not isinstance(Ts, list):
        Ts = [Ts]
    if not isinstance(Fs, list):
        Fs = [Fs]
    nx = get_backend(*Fs, *Ts, p)

    p = 1.0 / p
    tmpsum = sum(
        [
            lambdas[s] * nx.dot(Fs[s], Ts[s].T) * p[None, :]
            for s in range(len(Ts))
        ]
    )
    return tmpsum


def update_edge_feature_tensor(p, lambdas, T, Cs):
    r"""Updates :math:`\mathbf{C}` according to the l2 norm inner distance with
    the `S` :math:`\mathbf{T}_s` couplings.

    It is calculated at each iteration

    Parameters
    ----------
    p : array-like, shape (N,)
        Masses in the targeted barycenter.
    lambdas : list of float
        List of the `S` spaces' weights.
    T : list of S array-like of shape (ns,N)
        The `S` :math:`\mathbf{T}_s` couplings calculated at each iteration.
    Cs : list of S array-like, shape (ns,ns,d')
        Edge features tensors

    Returns
    -------
    C : array-like, shape (nt,nt,d')
        Updated :math:`\mathbf{C}` tensor.
    """
    p = list_to_array(p)
    T = list_to_array(*T)
    Cs = list_to_array(*Cs)
    if not isinstance(T, list):
        T = [T]
    if not isinstance(Cs, list):
        Cs = [Cs]
    nx = get_backend(*Cs, *T, p)

    # Proposition 2.10 in our paper
    tmpsum = sum(
        [
            lambdas[s]
            * nx.einsum(
                "ijd,jk...->ikd",
                nx.einsum("ij...,jkd->ikd", T[s].T, Cs[s]),
                T[s],
            )
            for s in range(len(T))
        ]
    )
    ppt = nx.reshape(nx.outer(p, p), shape=(len(p), len(p), 1))
    return tmpsum / ppt


def fused_network_gromov_wasserstein2(
    M,
    C1,
    C2,
    A1,
    A2,
    p,
    q,
    dist_fun_C="l2_norm",
    dist_fun_A="square_loss",
    alpha=0.33,
    beta=0.33,
    armijo=False,
    G0=None,
    log=False,
    **kwargs,
):
    r"""
    Computes the FNGW transport between two graphs

    See Algorithm 1 in our paper.

    Parameters
    ----------
    M : array-like, shape (ns, nt)
        Metric cost matrix between node features of source and target graphs
    C1 : array-like, shape (ns, ns, d')
        Edge feature tensor of the source graph
    C2 : array-like, shape (nt, nt, d')
        Edge feature tensor of the target graph
    A1 : array-like, shape (ns, ns)
        Structure matrix of the source graph
    A2 : array-like, shape (nt, nt)
        Structure matrix of the target graph
    p : array-like, shape (ns,)
        Distribution in the source space
    q : array-like, shape (nt,)
        Distribution in the target space
    dist_fun_C : str, optional
        Inner distance used for the edge feature tensor
    dist_fun_A : str, optional
        Loss function used for the structure matrix
    alpha : float, optional
        Trade-off parameter (0 < alpha < 1)
    beta : float, optional
        Trade-off parameter (0 < beta < 1)
    armijo : bool, optional
        If True the step of the line-search is found via an armijo research.
        Else closed form is used.
        If there are convergence issues use False.
    G0: array-like, shape (ns,nt), optional
        If None the initial transport plan of the solver is pq^T.
        Otherwise G0 must satisfy marginal constraints and will be used as
        initial transport of the solver.
    log : bool, optional
        record log if True
    **kwargs : dict
        parameters can be directly passed to the ot.optim.cg solver

    Returns
    -------
    fngw_dist : float
        FNGW distance for the given parameters.
    log : dict
        Log dictionary return only if log==True in parameters.
    """
    assert alpha + beta <= 1
    p, q = list_to_array(p, q)
    p0, q0, C10, C20, A10, A20, M0 = p, q, C1, C2, A1, A2, M
    if G0 is None:
        nx = get_backend(p0, q0, C10, C20, A10, A20, M0)
    else:
        G0_ = G0
        nx = get_backend(p0, q0, C10, C20, A10, A20, M0, G0_)

    p = nx.to_numpy(p)
    q = nx.to_numpy(q)
    C1 = nx.to_numpy(C10)
    C2 = nx.to_numpy(C20)
    A1 = nx.to_numpy(A10)
    A2 = nx.to_numpy(A20)
    M = nx.to_numpy(M0)

    if G0 is None:
        G0 = p[:, None] * q[None, :]
    else:
        G0 = nx.to_numpy(G0_)
        # Check marginals of G0
        np.testing.assert_allclose(G0.sum(axis=1), p, atol=1e-08)
        np.testing.assert_allclose(G0.sum(axis=0), q, atol=1e-08)

    constA, hA1, hA2 = init_matrix_A(A1, A2, p, q, dist_fun_A)
    constC = init_matrix_C(C1, C2, p, q, dist_fun_C)

    def f(G):
        return ngwloss(constC, C1, C2, G)

    def df(G):
        return ngwgrad(constC, C1, C2, G)

    def g(G):
        return gwloss(constA, hA1, hA2, G)

    def dg(G):
        return gwggrad(constA, hA1, hA2, G)

    T, cg_log = cg(
        p,
        q,
        (1 - alpha - beta) * M,
        reg_f=alpha,
        reg_g=beta,
        f=f,
        df=df,
        g=g,
        dg=dg,
        G0=G0,
        armijo=armijo,
        C1=C1,
        C2=C2,
        A1=A1,
        A2=A2,
        constC=constC,
        constA=constA,
        log=True,
        **kwargs,
    )

    fngw_dist = nx.from_numpy(cg_log["loss"][-1], type_as=C10)
    T0 = nx.from_numpy(T, type_as=C10)
    cg_log["fngw_dist"] = fngw_dist
    cg_log["u"] = nx.from_numpy(cg_log["u"], type_as=C10)
    cg_log["v"] = nx.from_numpy(cg_log["v"], type_as=C10)
    cg_log["T"] = T0

    # TODO: implement the gradient for p0, q0
    if dist_fun_C == "l2_norm" and dist_fun_A == "square_loss":
        gC1 = 2 * C1 * (p[:, None] * p[None, :])[:, :, None] - 2 * np.einsum(
            "ilt, kl->ikt", np.einsum("ij,jlt->ilt", T, C2), T
        )
        gC2 = 2 * C2 * (q[:, None] * q[None, :])[:, :, None] - 2 * np.einsum(
            "jkt, kl->jlt", np.einsum("ij,ikt->jkt", T, C1), T
        )
        gC1 = nx.from_numpy(gC1, type_as=C10)
        gC2 = nx.from_numpy(gC2, type_as=C10)

        gA1 = 2 * A1 * (p[:, None] * p[None, :]) - 2 * T.dot(A2).dot(T.T)
        gA2 = 2 * A2 * (q[:, None] * q[None, :]) - 2 * T.T.dot(A1).dot(T)
        gA1 = nx.from_numpy(gA1, type_as=A10)
        gA2 = nx.from_numpy(gA2, type_as=A10)

        fngw_dist = nx.set_gradients(
            fngw_dist,
            (p0, q0, C10, C20, A10, A20, M0),
            (
                cg_log["u"]
                - nx.mean(
                    cg_log["u"]
                ),  # No need for p0, q0 since they will not be updated,
                # keeps it right now
                cg_log["v"] - nx.mean(cg_log["v"]),
                alpha * gC1,
                alpha * gC2,
                beta * gA1,
                beta * gA2,
                (1 - alpha - beta) * T0,
            ),
        )
    if log:
        return fngw_dist, cg_log
    else:
        return fngw_dist


def init_matrix_C(C1, C2, p, q, dist="l2_norm"):
    r"""Computation of the sum of the first two terms of Equation (6) in our
    paper.

    Parameters
    ----------
    C1 : array-like, shape (ns, ns, d')
        Edge feature tensor of the source graph
    C2 : array-like, shape (nt, nt, d')
        Edge feature tensor of the target graph
    T :  array-like, shape (ns, nt)
        Coupling between source and target spaces
    p : array-like, shape (ns,)

    Returns
    -------
    constC : array-like, shape (ns, nt)

    """
    C1, C2, p, q = list_to_array(C1, C2, p, q)
    nx = get_backend(C1, C2, p, q)

    if dist == "l2_norm":

        def f1(a):
            return nx.sum(nx.power(a, 2), axis=-1)

        def f2(b):
            return nx.sum(nx.power(b, 2), axis=-1)

    else:
        raise ValueError

    constC1 = nx.dot(
        nx.dot(f1(C1), nx.reshape(p, (-1, 1))), nx.ones((1, len(q)), type_as=q)
    )
    constC2 = nx.dot(
        nx.ones((len(p), 1), type_as=p),
        nx.dot(nx.reshape(q, (1, -1)), f2(C2).T),
    )
    constC = constC1 + constC2

    return constC


def tensor_product(constC, C1, C2, T):
    r"""Implementation of the Prop. 2.5 in our paper.

    Parameters
    ----------
    constC : array-like, shape (ns, nt)
        the sum of the first two terms of Eq. (6)
    C1 : array-like, shape (ns, ns, d')
        Edge feature tensor of the source graph
    C2 : array-like, shape (nt, nt, d')
        Edge feature tensor of the target graph

    T : array-like, shape (ns, nt)

    Returns
    -------
    tens : array-like, shape (ns, nt)

    """
    constC, C1, C2, T = list_to_array(constC, C1, C2, T)
    nx = get_backend(constC, C1, C2, T)

    A = -2 * nx.einsum(
        "ijd, kjd->ikd", nx.einsum("ijd,jk...->ikd", C1, T), C2
    )  # (ns, nt, d)

    A = nx.sum(A, axis=-1)  # (ns, nt)
    tens = constC + A
    # tens -= tens.min()
    return tens


def ngwloss(constC, C1, C2, T):
    r"""Compute the third term of Eq.5 in our paper

    Parameters
    ----------
    constC : array-like, shape (ns, nt)
        the sum of the first two terms of Eq. (6)
    C1 : array-like, shape (ns, ns, d')
        Edge feature tensor of the source graph
    C2 : array-like, shape (nt, nt, d')
        Edge feature tensor of the target graph
    T : array-like, shape (ns, nt)
        Current value of transport matrix :math:`\mathbf{T}`
        Current value of transport matrix :math:`\mathbf{T}`

    Returns
    -------
    loss : float

    """

    tens = tensor_product(constC, C1, C2, T)

    tens, T = list_to_array(tens, T)
    nx = get_backend(tens, T)

    return nx.sum(tens * T)


def ngwgrad(constC, C1, C2, T):
    r"""Compute the third term of Eq.7 in our paper

    Parameters
    ----------
    constC : array-like, shape (ns, nt)
        the sum of the first two terms of Eq. (6)
    C1 : array-like, shape (ns, ns, d')
        Edge feature tensor of the source graph
    C2 : array-like, shape (nt, nt, d')
        Edge feature tensor of the target graph
    T : array-like, shape (ns, nt)
        Current value of transport matrix :math:`\mathbf{T}`

    Returns
    -------
    grad : array-like, shape (`ns`, `nt`)

    """
    return 2 * tensor_product(constC, C1, C2, T)


def cg(
    a,
    b,
    M,
    reg_f,
    reg_g,
    f,
    df,
    g,
    dg,
    G0=None,
    numItermax=200,
    numItermaxEmd=100000,
    stopThr=1e-9,
    stopThr2=1e-9,
    verbose=False,
    log=False,
    **kwargs,
):
    r"""
    Solve the general regularized OT problem with conditional gradient

        The function solves the following optimization problem:

    .. math::
        \gamma = \mathop{\arg \min}_\gamma \quad \langle \gamma, \mathbf{M}
        \rangle_F + \mathrm{reg} \cdot f(\gamma)

        s.t. \ \gamma \mathbf{1} &= \mathbf{a}

             \gamma^T \mathbf{1} &= \mathbf{b}

             \gamma &\geq 0
    where :

    - :math:`\mathbf{M}` is the (`ns`, `nt`) metric cost matrix
    - :math:`f` is the regularization term (and `df` is its gradient)
    - :math:`\mathbf{a}` and :math:`\mathbf{b}` are source and target weights
    (sum to 1)

    The algorithm used for solving the problem is conditional gradient as
    discussed in :ref:`[1] <references-cg>`


    Parameters
    ----------
    a : array-like, shape (ns,)
        samples weights in the source domain
    b : array-like, shape (nt,)
        samples in the target domain
    M : array-like, shape (ns, nt)
        loss matrix
    reg_f : float
        Regularization term >0
    reg_g : float
        Regularization term >0
    G0 :  array-like, shape (ns,nt), optional
        initial guess (default is indep joint density)
    numItermax : int, optional
        Max number of iterations
    numItermaxEmd : int, optional
        Max number of iterations for emd
    stopThr : float, optional
        Stop threshold on the relative variation (>0)
    stopThr2 : float, optional
        Stop threshold on the absolute variation (>0)
    verbose : bool, optional
        Print information along iterations
    log : bool, optional
        record log if True
    **kwargs : dict
             Parameters for linesearch

    Returns
    -------
    gamma : (ns x nt) ndarray
        Optimal transportation matrix for the given parameters
    log : dict
        log dictionary return only if log==True in parameters


    .. _references-cg:
    References
    ----------

    .. [1] Ferradans, S., Papadakis, N., Peyré, G., & Aujol, J. F. (2014).
    Regularized discrete optimal transport. SIAM Journal on Imaging Sciences,
    7(3), 1853-1882.

    See Also
    --------
    ot.lp.emd : Unregularized optimal ransport
    ot.bregman.sinkhorn : Entropic regularized optimal transport

    """
    a, b, M, G0 = list_to_array(a, b, M, G0)
    if isinstance(M, int) or isinstance(M, float):
        nx = get_backend(a, b)
    else:
        nx = get_backend(a, b, M)

    loop = 1

    if log:
        log = {"loss": []}

    if G0 is None:
        G = nx.outer(a, b)
    else:
        G = G0

    def cost(G):
        return nx.sum(M * G) + reg_f * f(G) + reg_g * g(G)

    f_val = cost(G)
    if log:
        log["loss"].append(f_val)

    it = 0

    if verbose:
        print(
            "{:5s}|{:12s}|{:8s}|{:8s}".format(
                "It.", "Loss", "Relative loss", "Absolute loss"
            )
            + "\n"
            + "-" * 48
        )
        print("{:5d}|{:8e}|{:8e}|{:8e}".format(it, f_val, 0, 0))

    while loop:

        it += 1
        old_fval = f_val

        # problem linearization
        Mi = M + reg_f * df(G) + reg_g * dg(G)
        # set M positive
        Mi += nx.min(Mi)

        # solve linear program
        Gc, logemd = emd(a, b, Mi, numItermax=numItermaxEmd, log=True)

        deltaG = Gc - G

        # line search
        alpha, fc, f_val = solve_linesearch(
            cost,
            G,
            deltaG,
            Mi,
            f_val,
            reg_f=reg_f,
            reg_g=reg_g,
            M=M,
            Gc=Gc,
            alpha_min=0.0,
            alpha_max=1.0,
            **kwargs,
        )

        G = G + alpha * deltaG

        # test convergence
        if it >= numItermax:
            loop = 0

        abs_delta_fval = abs(f_val - old_fval)
        relative_delta_fval = abs_delta_fval / abs(f_val)
        if relative_delta_fval < stopThr or abs_delta_fval < stopThr2:
            loop = 0

        if log:
            log["loss"].append(f_val)

        if verbose:
            if it % 20 == 0:
                print(
                    "{:5s}|{:12s}|{:8s}|{:8s}".format(
                        "It.", "Loss", "Relative loss", "Absolute loss"
                    )
                    + "\n"
                    + "-" * 48
                )
            print(
                "{:5d}|{:8e}|{:8e}|{:8e}".format(
                    it, f_val, relative_delta_fval, abs_delta_fval
                )
            )

    if log:
        log.update(logemd)
        return G, log
    else:
        return G


def solve_linesearch(
    cost,
    G,
    deltaG,
    Mi,
    val,
    armijo=True,
    C1=None,
    C2=None,
    reg_f=None,
    A1=None,
    A2=None,
    reg_g=None,
    Gc=None,
    constC=None,
    constA=None,
    M=None,
    alpha_min=None,
    alpha_max=None,
):
    """
    Solve the linesearch in the FNGW iterations

    Parameters
    ----------
    cost : method
        Cost in the FNGW for the linesearch
    G : array-like, shape(ns,nt)
        The transport map at a given iteration of the FNGW
    deltaG : array-like (ns,nt)
        Difference between the optimal map found by linearization in the FW
        algorithm and the value at a given iteration
    Mi : array-like (ns,nt)
        Cost matrix of the linearized transport problem. Corresponds to the
        gradient of the cost
    val : float
        Value of the cost at `G`
    armijo : bool, optional
        If True the steps of the line-search is found via an armijo research.
        Else closed form is used.
        If there is convergence issues use False.
    C1 : array-like (ns,ns,d'), optional
        Edge feature tensor of the source graph. Only used and necessary when
        armijo=False
    C2 : array-like (nt,nt,d'), optional
        Edge feature tensor in the target graph. Only used and necessary when
        armijo=False
    reg_f : float, optional
        Regularization parameter. Only used and necessary when armijo=False
    A1 : array-like (ns,ns), optional
        Structure matrix of the source graph. Only used and necessary when
        armijo=False
    A2 : array-like (nt,nt), optional
        Structure matrix of the target graph. Only used and necessary when
        armijo=False
    reg_g : float, optional
        Regularization parameter. Only used and necessary when armijo=False
    Gc : array-like (ns,nt)
        Optimal map found by linearization in the FW algorithm. Only used and
        necessary when armijo=False
    constC : array-like (ns,nt)
        Constant for the gromov cost. Only used and necessary when armijo=False
        See :ref:`[24] <references-solve-linesearch>`.

    M : array-like (ns,nt), optional
        Cost matrix between the features. Only used and necessary when
        armijo=False
    alpha_min : float, optional
        Minimum value for alpha
    alpha_max : float, optional
        Maximum value for alpha

    Returns
    -------
    alpha : float
        The optimal step size of the FW
    fc : int
        nb of function call. Useless here
    f_val : float
        The value of the cost for the next iteration

    """
    if armijo:
        # TODO: Update for armijo
        alpha, fc, f_val = line_search_armijo(
            cost, G, deltaG, Mi, val, alpha_min=alpha_min, alpha_max=alpha_max
        )
    else:
        G, deltaG, C1, C2, constC, A1, A2, constA, M = list_to_array(
            G, deltaG, C1, C2, constC, A1, A2, constA, M
        )
        if isinstance(M, int) or isinstance(M, float):
            nx = get_backend(G, deltaG, C1, C2, constC, A1, A2, constA)
        else:
            nx = get_backend(G, deltaG, C1, C2, constC, A1, A2, constA, M)

        dotC_1 = nx.sum(
            nx.einsum(
                "ijd, kjd->ikd", nx.einsum("ijd,jk...->ikd", C1, deltaG), C2
            ),
            axis=-1,
        )
        dotC_2 = nx.sum(
            nx.einsum("ijd, kjd->ikd", nx.einsum("ijd,jk...->ikd", C1, G), C2),
            axis=-1,
        )

        dotA_1 = nx.dot(nx.dot(A1, deltaG), A2.T)
        dotA_2 = nx.dot(nx.dot(A1, G), A2.T)

        a = -2 * reg_f * nx.sum(dotC_1 * deltaG) - 2 * reg_g * nx.sum(
            dotA_1 * deltaG
        )

        b = (
            nx.sum((M + reg_f * constC + reg_g * constA) * deltaG)
            - 2 * reg_f * (nx.sum(dotC_1 * G) + nx.sum(dotC_2 * deltaG))
            - 2 * reg_g * (nx.sum(dotA_1 * G) + nx.sum(dotA_2 * deltaG))
        )

        # c = cost(G)
        # c was pased to solve_linesearch as c, which does not exist
        alpha = solve_1d_linesearch_quad(a, b)
        if alpha_min is not None or alpha_max is not None:
            alpha = np.clip(alpha, alpha_min, alpha_max)
        fc = None
        f_val = cost(G + alpha * deltaG)

    return alpha, fc, f_val
