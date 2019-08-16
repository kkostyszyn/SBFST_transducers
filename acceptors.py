
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
                r = r_lt1
            elif f =="lt2":
                r = r_lt2
            elif f == "lt3":
                r = r_lt3
            else:
                r = repair
            x = (pynini.compose(s, r)).stringify()
            adv_list.append(x)         
        except:
            x = False
        
        t.write(s +  " | " + str(x) + " \n")

def alpha_gen(alph, e):
    zero = (e - e).optimize()
    sigma = zero
    for x in list(alph):
        sigma = A(x) | sigma
    sigma = pynini.closure(sigma)
    return sigma.optimize()

A = functools.partial(pynini.acceptor)
T = functools.partial(pynini.transducer)
e = pynini.epsilon_machine()

alpha = "abcd"
sigma = alpha_gen("abcd", e)
not_a = alpha_gen("bcd", e)
not_b = alpha_gen("acd", e)
not_ab = alpha_gen("cd", e)
sigma = sigma.optimize()

b = A("b")
a = A("a")

#------------------

repair = (not_b.star + T("b", "a") + sigma.star).optimize()
#for use with pt3
repair.write("r.fst")
repair_a = (not_a.star + T("a", "b") + sigma.star).optimize()
repair_a.write("r_a.fst")

###for use with lt1
### FIX LT1 TO MATCH NEW UNION
r_aaaa = (pynini.cdrewrite(T("a", "b"),
            "aaa", 
            "",
            sigma.star)).optimize()
r_bbbb = (pynini.cdrewrite(T("b", "a"),
            "bbb",
            "",
            sigma.star)).optimize()
r_lt1 = (pynini.compose(r_aaaa, r_bbbb)).optimize()
r_lt1.write("r_lt1.fst")
            
###for use with lt2
r_lt2 = (pynini.cdrewrite(T("bbbb", "bbba"),
                "",
                "",
                sigma.star)).optimize()
r_lt2.write("r_lt2.fst")

###for use with lt3
#delete one a - b^8, a^7
lt3_ab = pynini.cdrewrite(T("aaaaaaaa", "aaaaaaa"),
                        (not_a | "[BOS]"),
                        (not_a | "[EOS]"),
                        sigma)
#add one b - b^8 and a^7 
lt3_ba = pynini.cdrewrite(T("bbbbbbb", "bbbbbbbb"),
                        (not_b | "[BOS]"),
                        (not_b | "[EOS]"),
                        sigma)
r_lt3 = pynini.compose(lt3_ab, lt3_ba).optimize()
r_lt3.write("r_lt3.fst")

#------------------

#lt0 - bb
lt0_accept = (not_b.star + b + not_b.star + b + not_b.star).optimize()
lt0_accept.write("lt0_accept.fsa")

#lt1 - b^4 OR a^4
bbbb = (not_ab.star + b + b + b + b + not_ab.star).optimize()
aaaa = (not_ab.star + a + a + a + a + not_ab.star).optimize()
lt1_accept = pynini.union(aaaa, bbbb)
lt1_accept.write("lt1_accept.fsa")

#lt2 - b^4 AND a^4
aaaa_2 = (not_a.star + a + a + a + a + not_a.star).optimize()
bbbb_2 = (not_b.star + b + b + b + b + not_b.star).optimize()
lt2_accept = (pynini.intersect(aaaa_2, bbbb_2)).optimize()
lt2_accept.write("lt2_accept.fsa")

#lt3 - if b^8 then a^8 - NEEDS WORK 
#FIRST - strings with b^8 and a^8

b_8 = (not_b.star + b + b + b + b + b + b + b + b
            + not_b.star).optimize()
a_8 = (not_a.star + a + a + a + a + a + a + a + a
            + not_a.star).optimize()
b8_a8 = (pynini.intersect(b_8, a_8)).optimize()

#SECOND - strings with b^7 and a^7
b_7 = (not_b.star + b + b + b + b + b + b + b 
            + not_b.star).optimize()
a_7 = (not_a.star + a + a + a + a + a + a + a
            + not_a.star).optimize()
b7_a7 = (pynini.intersect(b_7, a_7)).optimize()
                        

lt3_accept = (pynini.union(b8_a8, b7_a7)).optimize()

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
