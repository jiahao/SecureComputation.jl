using SecureComputation, Test
SD = GF{982451653}
secret1 = 12345
secret2 = 54321

n = 16
t = 6

x = shamir(n, t, SD(secret1))
y = shamir(n, t, SD(secret2))

@test unshamir(x+y) == SD(secret1+secret2)
@test unshamir(x-x) == SD(0)
@test unshamir(x+(-x)) == SD(0)
@test unshamir(x*y) == SD(secret1*secret2)
@test unshamir(x/x) == SD(1)
@test unshamir(x\x) == SD(1)

SD = SecureComputation.GF{982451653}
n = 4
t = 6
x = rand(SD, n, n)
y = rand(SD, n, n)
z = x + y

sx = shamir.(t, t, x)
sy = shamir.(t, t, y)
sz = sx + sy

@test unshamir.(sz) == z

SD = GF{982451653}
n = 1
t = 3
x = rand(SD, n, n)
y = rand(SD, n, n)
z = x * y

sx = shamir.(2t+1, t, x)
sy = shamir.(2t+1, t, y)
sz = sx * sy

@test unshamir.(sz) == z
