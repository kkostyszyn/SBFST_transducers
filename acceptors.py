
import pynini
import functools

def list_string_set(ac):
    my_list = []
    paths = ac.paths()
    for s in paths.ostrings():
        my_list.append(s)
    my_list.sort(key=len)
    return my_list
    
def gen(fsa, accept):
    R = functools.partial(pynini.randgen)
    loop = 10
    n = 90
    for i in range(loop):
        num = int(n + n*i*0.1)
        rand = R(pynini.intersect(fsa, accept), npath=num, seed=0, select="uniform", max_length=1000, weighted=False)
    return list_string_set(rand)
    
def test(f, fsa, accept):
    
    t = open("results/" + f +"_results.txt", "w")
    
    adv_list = []
    for s in gen(fsa, accept):
        try:
            #change below to repair_a if testing pt3
            x = (pynini.compose(s, repair)).stringify()
            adv_list.append(x)         
        except:
            x = False
        
        t.write(s +  " | " + str(x) + " \n")
    



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

repair = (not_b.star + T("b", "a") + sigma.star).optimize()
#for use with pt3
repair_a = (not_a.star + T("a", "b") + sigma.star).optimize()


#------------------

#lt0 - bb
lt0_accept = (not_b.star + b + not_b.star + b + not_b.star).optimize()
lt0_accept.write("lt0_accept.fsa")

#lt1 - bbbb
lt1_accept = (not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star).optimize()
lt1_accept.write("lt1_accept.fsa")

#lt2, 5 * bb 
lt2_accept = (not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star).optimize()
lt2_accept.write("lt2_accept.fsa")

#lt3, 5 * bbbbbbbb
eight_bs = (not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star).optimize()
lt3_accept = (eight_bs + eight_bs + eight_bs + eight_bs + eight_bs).optimize()
lt3_accept.write("lt3_accept.fsa")


#pt0 - bb
pt0_accept = lt0_accept
pt0_accept.write("pt0_accept.fsa")


#pt1 - bbbbaaaa

pt1_accept = (not_b.star + b + 
                not_b.star + b + 
                not_b.star + b + 
                not_b.star + b + 
                not_b.star + a + 
                not_b.star + a + 
                not_b.star + a + 
                not_b.star + a + not_b.star).optimize()
pt1_accept.write("pt1_accept.fsa")


#pt2
pt2_accept = (not_b.star + b + 
                not_b.star + b + 
                not_b.star + b + 
                not_b.star + b + not_b.star).optimize()
pt2_accept.write("pt2_accept.fsa")


#pt3
pt3_accept = (not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + not_a.star).optimize()
pt3_accept.write("pt3_accept.fsa")



