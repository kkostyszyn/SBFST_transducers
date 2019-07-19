#
# A script to write out some simple FSAs
#

import pynini
import functools

A = functools.partial(pynini.acceptor) 
e = pynini.epsilon_machine()
zero = e - e
zero.optimize()

# Defining sigma and sigmastar
sigma4 = zero
for x in list("abcd"): sigma4 = A(x) | sigma4
sigma4.optimize()

sigma4Star = (sigma4.star).optimize()

a = A("a")
b = A("b")
c = A("c")
d = A("d")

# define simpler sigma and sigmastar
sigma2 = zero
for x in "ab": sigma2 = A(x) | sigma2
sigma2.optimize()
sigma2Star = sigma2.star.optimize()


def lg_containing_str(x,i):
    # return (sigma4Star + pynini.closure(b,i,i) + sigma4Star).minimize()
    return (sigma4Star + pynini.closure(x,i,i) + sigma4Star).minimize()


def lg_containing_ssq(x,i):
    return (pynini.closure(sigma4Star + x + sigma4Star,i,i)).minimize()


def lg_with_str(sigmastar, x, i):
    return (sigmastar + pynini.closure(x, i, i) + sigmastar).minimize()


def lg_with_ssq(sigmastar, x, i):
    return (pynini.closure(sigmastar + x + sigmastar, i, i)).minimize()


###############
# SL Examples #
###############

sl=dict()

sl[0] = sigma4Star - lg_containing_str(b,2)  # SL2 , forbidden bb
sl[1] = sigma4Star - lg_containing_str(b,4)  # SL4 , forbidden bbbb
sl[2] = sigma4Star - lg_containing_str(b,8)  # SL8 , forbidden bbbbbbbb

###############
# SP Examples #
###############

sp=dict()
sp[0] = sigma4Star - lg_containing_ssq(b,2)     # SP2 , forbidden bb
sp[1] = sigma4Star - lg_containing_ssq(b,4)     # SP4 , forbidden bbbb
sp[2] = sigma4Star - lg_containing_ssq(b,8)     # SP8 , forbidden bbbbbbbb

###############
# LT Examples #
###############

lt=dict()
# LT2 , at least one bb
lt[0] = lg_containing_str(b,2)

# LT4 , at least one bbbb or at least one aaaa
lt[1] = pynini.union(lg_containing_str(b,4), lg_containing_str(a,4))
# lt[1] = lg_containing_str(b,4) + lg_containing_str(a,4)

# LT4 , at least one bbbb and at least one aaaa
lt[2] = pynini.intersect(lg_containing_str(b,4), lg_containing_str(a,4))
        
# LT8 , if b^8 then a^8 (~~~ not b^8 or a^8)
lt[3] = (sigma4Star - lg_containing_str(b,8)) | lg_containing_str(a,8)

# aa and ab substrings
lt[4] = pynini.intersect(lg_containing_str(a, 2), lg_containing_str(a + b, 1))

# aa and ab substrings (using sigma = {a,b})
lt[5] = pynini.intersect(lg_with_str(sigma2Star, a, 2), lg_with_str(sigma2Star, a + b, 1))

###############
# PT Examples #
###############

pt=dict()
# PT2 , at least one bb
pt[0] = lg_containing_ssq(b,2)

# PT4 , at least one bbbb or at least one aaaa
pt[1] = pynini.union(lg_containing_ssq(b,4), lg_containing_ssq(a,4))
# pt[1] = lg_containing_ssq(b,4) + lg_containing_ssq(a,4)

# PT4 , at least one bbbb and at least one aaaa
pt[2] = pynini.intersect(lg_containing_ssq(b,4), lg_containing_ssq(a,4))
        
# PT8 , if b^8 then a^8 (~~~ not b^8 or a^8)
pt[3] = (sigma4Star - lg_containing_ssq(b,8)) | lg_containing_ssq(a,8)

# aa and ab subsequences
pt[4] = pynini.intersect(lg_containing_ssq(a, 2), sigma4Star + a + sigma4Star + b + sigma4Star)

# aa and ab subsequences (using sigma = {a,b})
pt[5] = pynini.intersect(lg_with_ssq(sigma2Star, a, 2), sigma2Star + a + sigma2Star + b + sigma2Star)


################
# LTT Examples #
################
ltt=dict()

# LTT t=2, k=1 "exactly two bs"
atleast2bs = (lg_containing_str(b,1) + lg_containing_str(b,1)).optimize()
forbid3bs = (sigma4Star - lg_containing_ssq(b,3)).optimize()
ltt[0] = atleast2bs * forbid3bs

# LTT t=2, k=2 "exactly two bb substrings"
atleast2bbs = (lt[0] + lt[0]).optimize()
forbid3bbs = sigma4Star - (atleast2bbs + lt[0]).optimize()
ltt[1] = atleast2bbs * forbid3bbs

# LTT t=5, k=2 "exactly five bb substrings"
atleast5bbs = pynini.closure(lt[0], 5, 5).optimize()
forbid6bbs = sigma4Star - pynini.closure(lt[0], 6, 6).optimize()
ltt[2] = atleast5bbs * forbid6bbs

# LTT t=5, k=8 "exactly five b^8 substrings"
atleast1b8 = pynini.closure(lg_containing_str(b,8), 5, 5).optimize()
forbid6b8s = sigma4Star - pynini.closure(lg_containing_str(b,8), 6, 6).optimize()
ltt[3] = atleast1b8 * forbid6b8s

ltt[2].write("lt2.fsa")
ltt[3].write("lt3.fsa")

