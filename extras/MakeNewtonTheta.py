import sympy 
import sys
from sympy.utilities.codegen import codegen
from sympy.printing import print_ccode

#The Data
M=sympy.Symbol('M')
m=sympy.Symbol('m')
E=sympy.Symbol('E')


#The Parameters
e=sympy.Symbol('a.error')
h=sympy.Symbol('a.freq')
F=sympy.Symbol('a.f')

#An array of the parameters.
params=[h,e,F]

#The three componenents of the likelihood function. H00 is homozygous for the major allele, H01 heterozygous, and H11 homozygous minor.

lnc=sympy.Symbol('lnc');
lne=sympy.Symbol('lne');
lneh=sympy.Symbol('lneh');
lnch=sympy.Symbol('lnch');
#H00=( h**2.+h*(1.-h)*F)*sympy.exp( lnc*M+lne*(m+e1+e2) ) 
#H01=2.*h*(1.-h)*(1.-F)*sympy.exp( lnch*(M+m)+lneh*(e1+e2) ) 
#H11=( (1.-h)**2.+h*(1.-h)*F)*sympy.exp( lnc*m+lne*(M+e1+e2) ) 

H00=( (1-h)**2.+h*(1.-h)*F)*( ( (1.-e)**M)*((e/3.)**(m+E) ) ) 
H01=2.*h*(1.-h)*(1.-F)*( ( ( (1.-e)/2.+(e/6.) )**(M+m) )*( (e/3.)**(E) ) )#-0.00001
H11=( (h)**2.+h*(1.-h)*F)*( ( (1.-e)**m )*( (e/3.)**(M+E) ) )#-0.00001 

#The log likelihood equation
lnL=sympy.log( H00+H01+H11 )

system_eq=[]

#We first need the three equations we are going to try and set to zero, i.e. the first partial derivitives wrt e h and F.
print "/*This code was automatically generated by "+str(sys.argv[0])+"*/\n"
print "#include \"allele_stat.h\""
print "#include \"quartet.h\""
print "#include \"typedef.h\""
print 
#numpy.set_printoptions(precission=18)
for x in range(0, 3):
	system_eq.append(sympy.diff(lnL, params[x]) )
	print "inline float_t H"+str(x)+" (const quartet_t &q, const allele_stat &a) {"
	print "\tconst float_t M=q.base[a.major];"
	print "\tconst float_t m=q.base[a.minor];"
	print "\tconst float_t E=q.base[a.e1]+q.base[a.e2];"
	sys.stdout.write("\treturn ")
	print_ccode(  system_eq[-1] )
	print ";\n}\n"

#Then we need to make the Jacobian, which is a matrix with ...
for x in range(0, 3):
	for y in range(0, 3):
		print "inline float_t J"+str(x)+str(y)+" (const quartet_t &q, const allele_stat &a) {"
		print "\tconst float_t M=q.base[a.major];"
		print "\tconst float_t m=q.base[a.minor];"
		print "\tconst float_t E=q.base[a.e1]+q.base[a.e2];"
		sys.stdout.write("\treturn ")
		print_ccode(sympy.diff(system_eq[x], params[y]) )
		print ";\n}\n"
print "inline float_t lnL_NR (const quartet_t &q, const allele_stat &a) {"
print "\tconst float_t M=q.base[a.major];"
print "\tconst float_t m=q.base[a.minor];"
print "\tconst float_t E=q.base[a.e1]+q.base[a.e2];"
sys.stdout.write("\treturn ")
print_ccode(sympy.simplify(lnL) )
print ";\n}\n"


quit()

