
import pynini
import functools


A = functools.partial(pynini.acceptor)
T = functools.partial(pynini.transducer)
e = pynini.epsilon_machine()

alpha = "abcd"
zero = (e - e).optimize()
sigma = zero
for x in list(alpha):
    sigma = A(x) | sigma
sigma = sigma.optimize()

b = A("b")
a = A("a")
not_b = (sigma - b | e).optimize()
not_a = (sigma - a | e).optimize()
sigma_star = (not_b | b).optimize()

#------------------

repair = (not_b.star + T("b", not_b) + sigma.star).optimize()

#------------------

#lt0 - bb
lt0_accept = (not_b.star + b + not_b.star + b + not_b.star).optimize()
lt0_accept.write("lt0_accept.fsa")

#lt1 - bbbb
lt1_accept = (not_b.star + b + not_b.star + b + not_b.star + b + not_b.star + b + not_b.star).optimize()
lt1_accept.write("lt1_accept.fsa")

#lt2, 3 seem identical to lt0? confirm 
lt2_accept = lt0_accept
lt2_accept.write("lt2_accept.fsa")

lt3_accept = lt0_accept
lt3_accept.write("lt3_accept.fsa")


#pt0 - bb
pt0_accept = lt0_accept
pt0_accept.write("pt0_accept.fsa")


#pt1 - bbbbaaaa

pt1_accept = (not_b.star + b + 
                not_b.star + b + 
                not_b.star + b + 
                not_b.star + b + 
                not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + pynini.intersect(not_b, not_a).star).optimize()
pt1_accept.write("pt1_accept.fsa")


#pt2
pt2_accept = (not_b.star + b + 
                not_b.star + b + 
                not_b.star + b + 
                not_b.star + b + not_b.star).optimize()
pt2_accept.write("pt2_accept.fsa")


#pt3
repair_a = (not_b.star + T("a", (sigma - a| e).optimize()) + sigma.star).optimize()
pt3_accept = (not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + not_a.star).optimize()
pt3_accept.write("pt3_accept.fsa")

