# find_primary_number.py
# Find the closest prime(s) to an integer entered on the console.
# Uses a fast deterministic Miller-Rabin witness set suitable for 64-bit integers,
# and falls back to probabilistic behavior for larger inputs (still very reliable).

def is_prime(n: int) -> bool:
    if n < 2:
        return False
    # small primes quick check
    small_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29)
    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False

    # Write n-1 as d*2^s1357267
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1

    # Deterministic Miller-Rabin bases for testing 64-bit integers
    witnesses = (2, 325, 9375, 28178, 450775, 9780504, 1795265022)

    def check(a: int) -> bool:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            return True
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                return True
        return False

    for a in witnesses:
        if a % n == 0:
            continue
        if not check(a):
            return False
    return True

def nearest_primes(n: int):
    if n >= 2 and is_prime(n):
        return [n], 0

    d = 0
    found_lower = None
    found_upper = None
    while True:
        d += 1
        lower = n - d
        upper = n + d

        if lower >= 2 and is_prime(lower):
            found_lower = lower
        if is_prime(upper):
            found_upper = upper

        if found_lower is not None or found_upper is not None:
            results = []
            if found_lower is not None:
                results.append(found_lower)
            if found_upper is not None:
                results.append(found_upper)
            return sorted(results), d

def prompt_int(prompt="Enter an integer: "):
    while True:
        try:
            return int(input(prompt))
        except ValueError:
            print("Please enter a valid integer.")

def main():
    n = prompt_int()
    primes, distance = nearest_primes(n)
    if distance == 0:
        print(f"{n} is prime.")
        print(f"Closest prime: {n} (distance 0)")
    else:
        if len(primes) == 1:
            print(f"Closest prime to {n}: {primes[0]} (distance {distance})")
        else:
            print(f"Two primes are equally close to {n}: {primes[0]} and {primes[1]} (distance {distance})")

if __name__ == "__main__":
    main()