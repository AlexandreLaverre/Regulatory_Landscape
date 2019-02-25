#!/usr/bin/env python3
# coding=utf-8



# Dict lifted interaction
inter_lift = {}

# Dict lifted_frag == real_frag
lift_real = {}

# Dict real interaction
inter_real = {}

# Conservation interaction
conserv_interaction = 0

for lifted_bait in inter_lift.keys():
    if lifted_bait in lift_real.keys():
        real_bait = lift_real[lifted_bait]
        for lifted_PIR in inter_lift[lifted_bait]:
            if lifted_PIR in lift_real.keys():
                real_PIR = lift_real[lifted_PIR]
                if real_PIR in inter_real[real_bait]:
                    conserv_interaction += 1


