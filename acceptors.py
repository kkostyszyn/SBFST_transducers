
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
    n = 100
    for i in range(loop):
        num = int(n + n*i*0.1)
        rand = R(pynini.intersect(fsa, accept), npath=num, seed=0, select="uniform", max_length=100, weighted=False)
    return list_string_set(rand)
    
def test(f, fsa, accept):
    #to run test - f must be the lang code, fsa must be the actual FSA in pynini (here found in the gen_FSAs/ folder, and accept must be the acceptor as defined below
    t = open("results/" + f +"_results.txt", "w")
    
    adv_list = []
    for s in gen(fsa, accept):
        try:
            #change below to repair_a if testing pt3
            if f == "pt3":
                r = repair_a
            elif f == "lt1":
                r = r_aaaa
            elif f =="lt2":
                r = repair_lt2
            else:
                r = repair
            x = (pynini.compose(s, r)).stringify()
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
sigma = pynini.closure(sigma)
sigma = sigma.optimize()

b = A("b")
a = A("a")
not_b = (sigma - b | e).optimize()
not_a = (sigma - a | e).optimize()

#------------------

repair = (not_b.star + T("b", "a") + sigma.star).optimize()
#for use with pt3
repair_a = (not_a.star + T("a", "b") + sigma.star).optimize()

###for use with lt1
###when using lt1_accept = aaaa, use r_aaaa. when using lt1_accept = bbbb, use r_bbbb
r_aaaa = (pynini.cdrewrite(T("a", "b"),
            "aaa", 
            "",
            sigma.star)).optimize()
r_bbbb = (pynini.cdrewrite(T("b", "a"),
            "bbb",
            "",
            sigma.star)).optimize()
            
###for use with lt2
repair_lt2 = (pynini.cdrewrite(T("bbbb", "bbba"),
                "",
                "",
                sigma.star)).optimize()


#------------------

#lt0 - bb
lt0_accept = (not_b.star + b + not_b.star + b + not_b.star).optimize()
lt0_accept.write("lt0_accept.fsa")

#lt1 - b^4 OR a^4
bbbb = (not_b.star + b + b + b + b + not_b.star).optimize()
aaaa = (not_a.star + a + a + a + a + not_a.star).optimize()
lt1_accept = aaaa
lt1_accept.write("lt1_accept.fsa")

#lt2 - b^4 AND a^4
#use aaaa and bbbb above
lt2_accept = (not_b.star + b + b + b + b + not_b.star).optimize()
lt2_accept.write("lt2_accept.fsa")

#lt3 - if b^8 then a^8
lt3_accept = (not_b.star + b + b + b + b + b + b + b + b
            + not_b.star).optimize()
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


#pt2 - bbbb
pt2_accept = (not_b.star + b + 
                not_b.star + b + 
                not_b.star + b + 
                not_b.star + b + not_b.star).optimize()
pt2_accept.write("pt2_accept.fsa")


#pt3 - aaaaaaaa
pt3_accept = (not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + 
                not_a.star + a + not_a.star).optimize()
pt3_accept.write("pt3_accept.fsa")

################
# IGNORE THESE #
################

#these were LTT mistakenly made as LT

#ltt1 - bbbb
ltt1_accept = (not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star).optimize()
ltt1_accept.write("ltt1_accept.fsa")

#ltt2, 5 * bb 
ltt2_accept = (not_b.star + b 
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
#ltt2_accept.write("ltt2_accept.fsa")

#ltt3, 5 * bbbbbbbb
eight_bs = (not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star + b 
            + not_b.star + b).optimize()
ltt3_accept = (eight_bs + eight_bs + eight_bs + eight_bs + eight_bs + not_b.star).optimize()
#ltt3_accept.write("ltt3_accept.fsa")
