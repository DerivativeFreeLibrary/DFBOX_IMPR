#============================================================================================
#    SDBOX - FORTRAN90 implementation of a Derivative-Free algorithm for bound 
#    constrained optimization problems 
#    Copyright (C) 2011  G.Liuzzi, S. Lucidi
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    S. Lucidi, M. Sciandrone. A Derivative-Free Algorithm for Bound Constrained Optimization, 
#    Computational Optimization and Applications, 21(2): 119-142 (2002)
#    DOI: 10.1023/A:1013735414984
#
#============================================================================================
import math
import time
import problem
funct = problem.funct
setbounds = problem.setbounds
startp = problem.startp

#import problem1
#funct = problem1.funct
#setbounds = problem1.setbounds
#startp = problem1.startp

sqrt = math.sqrt

#     *********************************************************
#     *         
#     *                 stop condition verify
#     *
#     *********************************************************
def stop(n,alfa_d,nf,ni,fstop,f,alfa_stop,nf_max, flag_fail):

	istop = 0
	alfa_max = 0.0

	for i in range(0,n):
		if alfa_d[i] > alfa_max:
			alfa_max = alfa_d[i]
     
	if ni >= (n+1):
		ffm = f
		for i in range(0,n):
			ffm += fstop[i]

		ffm /= float(n+1)

		ffstop = (f-ffm)**2

		for i in range(0,n):
			ffstop += (fstop[i]-ffm)**2

		ffstop = sqrt(ffstop/float(n+1))

      
	if alfa_max <= alfa_stop:
		istop = 1

	if nf > nf_max:
		istop = 2

	return istop, alfa_max

#     *********************************************************
#     *         
#     *                 Continuous Linesearch
#     *
#     *********************************************************
def linesearchbox_cont(n,x,f,d,alfa_d,j,alfa_max,iprint,bl,bu,nf):
      
	z = [a for a in x]
	  
	gamma = 1.e-6
	delta =0.5
	delta1=0.5
	i_corr_fall =0
	ifront      =0

	if iprint >= 1:
		print('j =%d    d(j) =%f alfa=%e' % (j,d[j-1],alfa_d[j-1]))

	if abs(alfa_d[j-1]) <= 1.e-3*min(1.0,alfa_max):
		alfa = 0.0
		if iprint >= 1:
			print('  alfa piccolo')
			print(' alfa_d(j)=%e    alfamax=%e' % (alfa_d[j-1],alfa_max))
		return alfa, f, nf, i_corr_fall
      
	for ielle in range(1,3):

		if d[j-1] > 0.0:

			if (alfa_d[j-1]-(bu[j-1]-x[j-1])) < -1.e-6:
				alfa = max(1.e-24,alfa_d[j-1])
			else:
				alfa = bu[j-1]-x[j-1]
				ifront=1
				if iprint >= 1:
					print(' point on the boundary. *')
		else:

			if (alfa_d[j-1]-(x[j-1]-bl[j-1])) < -1.e-6:
				alfa = max(1.e-24,alfa_d[j-1])
			else:
				alfa = x[j-1]-bl[j-1]
				ifront=1
				if iprint >= 1:
					print(' point on the boundary. *')

		if abs(alfa) <= 1.e-3*min(1.0,alfa_max):

			d[j-1] = -d[j-1]
			i_corr_fall += 1
			ifront = 0

			if iprint >= 1:
				print(' direzione opposta per alfa piccolo')
				print(' j =%d    d(j) =%f' % (j,d[j-1]))
				print(' alfa=%e    alfamax=%e' % (alfa,alfa_max))

			alfa = 0.0
			continue

		alfaex = alfa
		z[j-1] = x[j-1] + alfa*d[j-1]

		fz  = funct(z)
		nf += 1

		if iprint >= 1:
			print(' fz =%f   alfa =%e' % (fz,alfa))

		if iprint >= 2:
			for i in range(0,n):
				print(' z(%d)=%f' % (i,z[i]))

		fpar = f - gamma*alfa**2

		if fz < fpar:

		# expansion step

			while True:

				if ifront==1:

					if iprint >= 1:
						print(' accetta punto sulla frontiera fz =%f   alfa =%f' % (fz,alfa))

					alfa_d[j-1] = delta*alfa
					return alfa, fz, nf, i_corr_fall

				if d[j-1] > 0.0:

					if (alfa/delta1-(bu[j-1]-x[j-1])) < -1.e-6:
						alfaex = alfa/delta1
					else:
						alfaex = bu[j-1]-x[j-1]
						ifront = 1
						if iprint >= 1:
							print(' punto espan. sulla front.')

				else:

					if (alfa/delta1-(x[j-1]-bl[j-1])) < -1.e-6:
						alfaex = alfa/delta1
					else:
						alfaex = x[j-1]-bl[j-1]
						ifront = 1
						if iprint >= 1:
							print(' punto espan. sulla front.')

				z[j-1] = x[j-1] + alfaex*d[j-1] 

				fzdelta = funct(z)
				nf     += 1

				if iprint >= 1:
					print(' fzex=%f  alfaex=%f' % (fzdelta,alfaex))

				fpar = f - gamma*alfaex**2

				if fzdelta < fpar:
					fz   = fzdelta
					alfa = alfaex
				else:
					alfa_d[j-1] = delta*alfa

					if iprint>= 1:
						print(' accetta punto fz =%f   alfa =%f' % (fz,alfa))

					return alfa, fz, nf, i_corr_fall

		else:   #opposite direction    

			d[j-1] = -d[j-1]
			ifront = 0

			if iprint >= 1:
				print(' direzione opposta')
				print(' j =%d    d(j) =%f' % (j,d[j-1]))

	if i_corr_fall != 2:
		alfa_d[j-1] = delta*alfa_d[j-1]

	alfa = 0.0

	if iprint >= 1:
		print(' failure along the direction')

	return alfa, f, nf, i_corr_fall

#     *********************************************************
#     *         
#     *                 DF_BOX outer loop
#     *
#     *********************************************************
def sd_box(n,x,f,bl,bu,alfa_stop,nf_max,maxiter,nf,iprint):
#     initialization
	eta       = 1.e-6
	num_fal   = 0
	istop     = 0
	flag_fail = [0]*n
	fstop     = [0.0]*(n+1)
	alfa_d    = [0.0]*n
	d         = [1.0]*n
	
	format100 = ' ni=%4d  nf=%5d   f=%12.5e   alfamax=%12.5e'

#---- choice of the starting stepsizes along the directions --------

	for i in range(0,n):
		alfa_d[i] = max(1.e-3,min(1.0,abs(x[i])))

		if iprint >= 1:
			print(' alfainiz(%d)=%e' % (i,alfa_d[i]))

	alfa_max = max(alfa_d)
	f   = funct(x)
	nf += 1
	i_corr = 1
	fstop[i_corr-1] = f

#---------------------------   
#     main loop
#---------------------------

	for ni in range(1,maxiter+1):

		if iprint >= 0:
			print(format100 % (ni,nf,f,alfa_max))

#-------------------------------------
#    sampling along coordinate i_corr
#-------------------------------------
		alfa, fz, nf, i_corr_fall = linesearchbox_cont(n,x,f,d,alfa_d,i_corr,alfa_max,iprint,bl,bu,nf)

		if abs(alfa) >= 1.e-12:
			flag_fail[i_corr-1] = 0
			x[i_corr-1] = x[i_corr-1] + alfa*d[i_corr-1]
			f = fz
			fstop[i_corr-1] = f
			num_fal = 0
			ni += 1

		else:

			flag_fail[i_corr-1] = 1
			if i_corr_fall < 2:
				fstop[i_corr-1] = fz
				num_fal += 1
				ni += 1

		if i_corr < n:
			i_corr += 1
		else:
			i_corr = 1

		istop, alfa_max = stop(n,alfa_d,nf,ni,fstop,f,alfa_stop,nf_max,flag_fail)

		if istop >= 1: 
			if iprint >= 0:
				print(format100 % (ni,nf,f,alfa_max))
			break
	return x, f


#-----------------------------------------------------------------------
#      Starting point and bound calculation
#-----------------------------------------------------------------------

bl, bu = setbounds()	
x = startp()
n = len(x)

exit = False

for i in range(0,n):
	if (x[i] < bl[i]) or (x[i] > bu[i]):
		print('ERROR: initial point is out of the feasible box')
		exit = True		

if not exit:
	num_funct   = 0 
	alfa_stop   = 1.e-6
	maxiter     = 20000
	nf_max      = 20000
	iprint      = 0

	fob   = funct(x)
	finiz = fob
	num_funct += 1

	print(' ------------------------------------------------- ')
	print(' objective function at xo = %f' % fob)
	print(' ------------------------------------------------- ')

	tbegin = time.clock()

	xstar, f = sd_box(n,x,fob,bl,bu,alfa_stop,nf_max,maxiter,num_funct,iprint)

	tend = time.clock()

	format = ' & %3d & %14.7e & %14.7e & %5d & %9.2e'
	print(format % (n,finiz,f,num_funct,(tend-tbegin)))


	print('------------------------------------------------------------------------------')
	print(' total time:%f' % (tend-tbegin))
	print('------------------------------------------------------------------------------')
