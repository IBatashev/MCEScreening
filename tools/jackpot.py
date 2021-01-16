import random
import time

poker = ['V', 'o', '@', '*', '#', '!', '7']
print('\n')
speed = 100
for a in range(1,speed):
    if a < speed*0.6: i = random.randrange(7)
    else: i = 7
    if a < speed*0.8: j = random.randrange(7)
    else: j = 7
    if a < speed-1: k = random.randrange(7)
    else: k = 7
    time.sleep(0.07)
    print('\r', '\t\t[| {} || {} || {} |]'.format(*[poker[i-1], poker[j-1], poker[k-1]]), end='')
time.sleep(0.7)
print("\n\n\t\t     Jackpot!")
