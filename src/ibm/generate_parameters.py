#!/usr/bin/env python


exe_name = "xbadger"

tau = 2.0
phi = 1.0

mortality_rate = 0.01
mortality_rate_cubs_alone = 0.9
mortality_rate_cubs_highK = 0.9

K = 20;

# value at which virulence is 0.5 that of maximum
tau = 0.01

# transmission coefficient
phi = 1.0

mu = 0.25
sdmu = 0.01

init_v =0.5


nrep = 1

counter = 1

for i in range(0,nrep):


    print("echo " + str(counter))

    counter += 1 

    print("./" + exe_name + " " + \
            str(tau) + " " +\
            str(phi) + " " +\
            str(mortality_rate) + " " +\
            str(mortality_rate_cubs_alone) + " " +\
            str(mortality_rate_cubs_highK) + " " +\
            str(init_v) + " " +\
            str(K) + " " +\
            str(mu) + " " +\
            str(sdmu) + " "
            )

