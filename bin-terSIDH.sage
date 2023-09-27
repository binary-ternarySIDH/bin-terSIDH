import sys
import random
import time

from montgomery_xz import KummerLine
from montgomery_isogeny import KummerLineIsogenyComposite
from utilities import torsion_basis

proof.all(False)
random.seed(int(42))

lam = 128
ternary = True

for arg in sys.argv[1:]:
    if arg.lower() in ["--192"]:
        lam = 192
        break
    elif arg.lower() in ["--256"]:
        lam = 256
        break

for arg in sys.argv[1:]:
    if arg.lower() in ["--bin", "--binary"]:
        ternary = False
        break

params_binSIDH = [134, 192, 256]
params_terSIDH = [93, 128, 162]

if ternary:
    params = params_terSIDH
    sk_choices = [0, 1, 2]
else:
    params = params_binSIDH
    sk_choices = [1, 2]

if lam == 128:
    t = params[0]
elif lam == 192:
    t = params[1]
elif lam == 256:
    t = params[2]
else:
    raise Exception("The security parameter needs to be 128, 192, or 256.")




def make_prime(p):
    '''
    Given a value `p`, finds a cofactor `f`
    such that p*f - 1 is prime.
    '''
    for i in range(1000):
        if (p*i - 1).is_prime():
            return p*i - 1, i

def compute_kernel_scalars(s, Alice=True):
    """ 
    Given a ternary secret `s`, returns scalars `B0` and `B1`
    such that the isogeny associated with `s` and orientation (P, Q)
    has kernel <[B0]*P + [B1]*Q>.
    The function also returns `order0` and `order1`, the orders
    of points [B0]*P and [B1]*Q, which is used in the isogeny computations.
    """
    B0 = 1
    B1 = 1
    order0 = 1
    order1 = 1

    t = len(s)

    if Alice:
        P = Primes()[:2*t:2]
    else:
        P = Primes()[1:2*t:2]

    for i, p in enumerate(P):
        if Alice and i == 0:
            p = 4

        if s[i] == 2:
            B1 *= p
            order0 *= p
        elif s[i] == 1:
            B0 *= p
            order1 *= p
        else:
            B0 *= p
            B1 *= p
    
    return B0, B1, order0, order1


##### Setup ############################

A = 2*prod(Primes()[:2*t:2]) # The 2* ensures that p \equiv 3 mod 4
B = prod(Primes()[1:2*t:2])

p, f = make_prime(A*B)


FF.<x> = GF(p)[]
F.<i> = GF(p^2, modulus=x^2+1) 
E0 = EllipticCurve(F, [0, 6, 0, 1, 0])
E0.set_order((p+1)**2) 

PA, QA = torsion_basis(E0, A)
PB, QB = torsion_basis(E0, B)

## Ensures that 2*PA != (0,0) and
## 2*QA != (0,0), which causes problems in
## computing isogenies with x-only arithmetic
if A//2 * PA == E0(0, 0):
    PA = PA + QA
elif A//2 * QA == E0(0, 0):
    QA = PA + QA 

assert A//2 * PA != E0(0, 0) and A//2 * QA != E0(0, 0)


E0 = KummerLine(E0)
xPA, xQA = E0(PA[0]), E0(QA[0])
xPB, xQB = E0(PB[0]), E0(QB[0])



def keygen(A_points, B_points, Alice=True):
    # Extract Kummer points
    xP, xQ = A_points
    xR, xS = B_points

    # Generate the secret data
    sk = random.choices(sk_choices, k=t)
    B0, B1, order0, order1 = compute_kernel_scalars(sk, Alice=Alice)
    sk = (B0, B1, order0, order1)
    
    # Compute the isogeny kernels
    xK0 = xP * B0
    xK1 = xQ * B1

    # Compute the first isogeny
    phi0 = KummerLineIsogenyComposite(E0, xK0, order0)
    xK1 = phi0(xK1)

    # Evaluate action on aux data
    xR, xS = phi0(xR), phi0(xS)

    # Compute the second isogeny from the codomain of phi0
    E1 = phi0.codomain()
    phi1 = KummerLineIsogenyComposite(E1, xK1, order1)

    # Evaluate action on aux data
    xR, xS = phi1(xR), phi1(xS)

    # Generate the masking values
    modulus = B if Alice else A # modulus is the order of torsion points
    mask = modulus
    while gcd(mask, modulus) != 1:
        mask = random.randint(0, modulus)
    mask_inv = 1/mask % modulus

    # Scale the torsion images
    xR = mask * xR
    xS = mask_inv * xS

    # Compute the public data
    E1 = phi1.codomain()
    pk = (E1, xR, xS)

    return sk, pk


def shared(sk, pk):
    # Parse the private/public keys
    B0, B1, order0, order1 = skA
    E, xR, xS = pkB
    
    # Compute the isogeny kernels
    xK0 = xR * B0
    xK1 = xS * B1

    # Compute the first isogeny
    phi0 = KummerLineIsogenyComposite(E, xK0, order0)
    xK1 = phi0(xK1)

    # Compute the second isogeny from the codomain of phi0
    E = phi0.codomain()
    phi1 = KummerLineIsogenyComposite(E, xK1, order1)

    EAB = phi1.codomain()
    return EAB.j_invariant()


N = 1
tt = [0, 0, 0, 0]

for _ in range(N):
    t0 = time.process_time_ns()
    skA, pkA = keygen((xPA, xQA), (xPB, xQB), Alice=True)
    tt[0] += time.process_time_ns() - t0

    t0 = time.process_time_ns()
    skB, pkB = keygen((xPB, xQB), (xPA, xQA), Alice=False)
    tt[1] += time.process_time_ns() - t0

    t0 = time.process_time_ns()
    ssA = shared(skA, pkB)
    tt[2] += time.process_time_ns() - t0

    t0 = time.process_time_ns()
    ssB = shared(skB, pkA)
    tt[3] += time.process_time_ns() - t0

    assert ssA == ssB

tt = [float(t) / N / 10^6 for t in tt]

print(f"KeyGen_A took {(tt[0]):.1f} ms")
print(f"KeyGen_B took {(tt[1]):.1f} ms")
print(f"shared_A took {(tt[2]):.1f} ms")
print(f"shared_B took {(tt[3]):.1f} ms")