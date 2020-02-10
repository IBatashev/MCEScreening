import random
import time

poker = ['V', 'o', '@', '*', '#', '!', '7']
print('\n')
for a in range(1,201):
    if a < 120: i = random.randrange(7)
    else: i = 7
    if a < 160: j = random.randrange(7)
    else: j = 7
    if a < 200: k = random.randrange(7)
    else: k = 7
    time.sleep(0.07)
    print('\r', '\t\t[| {} || {} || {} |]'.format(*[poker[i-1], poker[j-1], poker[k-1]]), end='')
time.sleep(0.7)
print("\n\n\t\t     Jackpot!")