
import numpy as np
import scipy as sp
import math
import pandas as pd 
import time
import geweke
import os
import gc
import sys


def sample_gamma(beta,sigma_0,sigma_1,pie):
	p = np.empty(len(beta))
	d1 = pie*sp.stats.norm.pdf(beta,loc=0,scale=sigma_1)
	d0 = (1-pie)*sp.stats.norm.pdf(beta,loc=0,scale=sigma_0)
	p = d1/(d0+d1)
	gamma = np.random.binomial(1,p)
	return(gamma)

def sample_pie(gamma,pie_a,pie_b):
	a_new = np.sum(gamma)+pie_a
	b_new = np.sum(1-gamma)+pie_b
	pie_new = np.random.beta(a_new,b_new)
	return(pie_new)

def sample_sigma_1(beta,gamma,a_sigma,b_sigma):
	a_new = 0.5*np.sum(gamma)+a_sigma
	b_new = 0.5*np.sum(np.multiply(np.square(beta),gamma))+b_sigma
	sigma_1_neg2 =np.random.gamma(a_new,1.0/b_new)
	sigma_1_new = math.sqrt(1/sigma_1_neg2)
	return(sigma_1_new)

def sample_sigma_e(y,H_beta,C_alpha,a_e,b_e):
	n = len(y)
	a_new = float(n)/2+a_e
	resid = y - H_beta - C_alpha
	b_new = np.sum(np.square(resid))/2+b_e
	sigma_e_neg2 =np.random.gamma(a_new,1.0/b_new)
	sigma_e_new = math.sqrt(1/sigma_e_neg2)
	return(sigma_e_new)

def sample_alpha(y,H_beta,C_alpha,C,alpha,sigma_e,C_norm_2):

	r,c = C.shape

	if c == 1:
		#new_variance = 1/(np.linalg.norm(C[:,0])**2*sigma_e**-2)
		new_variance = 1/(C_norm_2[0]*sigma_e**-2)
		new_mean = new_variance*np.dot((y-H_beta),C[:,0])*sigma_e**-2
		alpha = np.random.normal(new_mean,math.sqrt(new_variance))
		C_alpha = C[:,0] * alpha
	else:
		for i in range(c):
			#new_variance = 1/(np.linalg.norm(C[:,i])**2*sigma_e**-2)
			new_variance = 1/(C_norm_2[i]*sigma_e**-2)
			C_alpha_negi = C_alpha - C[:,i] * alpha[i]
			new_mean = new_variance*np.dot(y-C_alpha_negi-H_beta,C[:,i])*sigma_e**-2
			alpha[i] = np.random.normal(new_mean,math.sqrt(new_variance))
			C_alpha = C_alpha_negi + C[:,i] * alpha[i]

	return(alpha,C_alpha)

def sample_beta(y,C_alpha,H_beta,H,beta,gamma,sigma_0,sigma_1,sigma_e,H_norm_2):

	sigma_e_neg2 = sigma_e**-2
	sigma_0_neg2 = sigma_0**-2
	sigma_1_neg2 = sigma_1**-2

	for i in range(len(beta)):
		H_beta_negi = H_beta - H[:,i] * beta[i]
		residual = y - C_alpha -  H_beta + H[:,i] * beta[i]
		new_variance = 1/(H_norm_2[i]*sigma_e_neg2+(1-gamma[i])*sigma_0_neg2+gamma[i]*sigma_1_neg2)
		new_mean = new_variance*np.dot(residual,H[:,i])*sigma_e_neg2
		beta[i] = np.random.normal(new_mean,math.sqrt(new_variance))
		H_beta = H_beta_negi + H[:,i] * beta[i]

	return(beta,H_beta)

def sampling(verbose,y,C,HapDM,sig0_initiate,iters,prefix,num,trace_container,gamma_container,beta_container,alpha_container):

	## set random seed for the process
	np.random.seed(int(time.time()) + os.getpid())

	#initiate beta,gamma and H matrix
	H = np.array(HapDM)

	H_r,H_c = H.shape


	C_c = C.shape[1]


	##specify hyper parameters
	pie_a = 1
	pie_b = H_c / 100
	a_sigma = 1
	b_sigma = 1
	a_e = 1
	b_e = 1

	sigma_0 = sig0_initiate
	sigma_1 = math.sqrt(1/np.random.gamma(a_sigma,b_sigma))
	sigma_e = math.sqrt(1/np.random.gamma(a_e,b_e))
	pie = np.random.beta(pie_a,pie_b)

	print("initiation for chain %i:" %(num) ,sigma_1,sigma_e,pie)
	
	#initiate alpha, alpha_trace, beta_trace and gamma_trace

	it = 0
	burn_in_iter = 2000
	trace = np.empty((iters-2000,5))
	alpha_trace = np.empty((iters-2000,C_c))
	gamma_trace = np.empty((iters-2000,H_c))
	beta_trace = np.empty((iters-2000,H_c))
	top5_beta_trace = np.empty((iters-2000,5))

	alpha = np.random.random(size = C_c)
	gamma = np.random.binomial(1,pie,H_c)
	beta = np.array(np.zeros(H_c))

	for i in range(H_c):
		if gamma[i] == 0:
			beta[i] = np.random.normal(0,sigma_0)
		else:
			beta[i] = np.random.normal(0,sigma_1) 

	#start sampling

	H_beta = np.matmul(H,beta)
	C_alpha = np.matmul(C,alpha)


	## precompute some variables 

	C_norm_2 = np.sum(C**2,axis=0)
	H_norm_2 = np.sum(H**2,axis=0)

	while it < iters:
		before = time.time()
		sigma_1 = sample_sigma_1(beta,gamma,a_sigma,b_sigma)
		if sigma_1 < 0.05:
			sigma_1 = 0.05
			pie = 0
		else:
			pie = sample_pie(gamma,pie_a,pie_b)
		sigma_e = sample_sigma_e(y,H_beta,C_alpha,a_e,b_e)
		gamma = sample_gamma(beta,sigma_0,sigma_1,pie)
		alpha,C_alpha = sample_alpha(y,H_beta,C_alpha,C,alpha,sigma_e,C_norm_2)
		beta,H_beta = sample_beta(y,C_alpha,H_beta,H,beta,gamma,sigma_0,sigma_1,sigma_e,H_norm_2)
		genetic_var = np.var(H_beta)
		pheno_var = np.var(y - C_alpha)
		large_beta = np.absolute(beta) > 0.3
		large_beta_ratio = np.sum(large_beta) / len(beta)
		total_heritability = genetic_var / pheno_var

		after = time.time()
		if it > 100 and total_heritability > 1:
			continue

		else:
			if verbose:
				print(num,it,str(after - before),sigma_1,sigma_e,large_beta_ratio,total_heritability,sum(gamma))

			if it >= burn_in_iter:
				trace[it-burn_in_iter,:] = [sigma_1,sigma_e,large_beta_ratio,total_heritability,sum(gamma)]
				gamma_trace[it-burn_in_iter,:] = gamma
				beta_trace[it-burn_in_iter,:] = beta
				alpha_trace[it-burn_in_iter,:] = alpha
				top5_beta_trace[it-burn_in_iter,:] = np.sort(np.absolute(beta))[::-1][:5]

			if it >= burn_in_iter + 9999: # after burn-in iterations, test convergence

				max_z = []
				
				for a in range(C_c):
					after_burnin_alpha = alpha_trace[:,a]
					alpha_zscores = geweke.geweke(after_burnin_alpha)[:,1]
					max_z.append(np.amax(np.absolute(alpha_zscores)))

				for b in range(5):
					after_burnin_beta = top5_beta_trace[:,b]
					beta_zscores = geweke.geweke(after_burnin_beta)[:,1]
					max_z.append(np.amax(np.absolute(beta_zscores)))

				#convergence for large beta ratio
				after_burnin_pie = trace[:,2]
				pie_zscores = geweke.geweke(after_burnin_pie)[:,1]
				max_z.append(np.amax(np.absolute(pie_zscores)))

				#convergence for total heritability
				after_burnin_var = trace[:,4]
				var_zscores = geweke.geweke(after_burnin_var)[:,1]
				max_z.append(np.amax(np.absolute(var_zscores)))

				# #convergence for sigma_1
				# after_burnin_sigma1 = trace[:,1]
				# sigma1_zscores = geweke.geweke(after_burnin_sigma1)[:,1]
				# max_z.append(np.amax(np.absolute(sigma1_zscores)))

				#convergence for sigma_e
				after_burnin_sigmae = trace[:,1]
				sigmae_zscores = geweke.geweke(after_burnin_sigmae)[:,1]
				max_z.append(np.amax(np.absolute(sigmae_zscores)))
				
				if  np.amax(max_z) < 1.5:
					print("Chain %i: convergence has been reached at %i iterations." %(num,it))
					break

				else:
					trace_ = np.empty((1000,5))
					gamma_trace_ = np.empty((1000,H_c))
					beta_trace_ = np.empty((1000,H_c))
					alpha_trace_ = np.empty((1000,C_c))
					top5_beta_trace_ = np.empty((1000,5))

					trace = np.concatenate((trace[-(iters - burn_in_iter-1000):,:],trace_),axis=0)
					gamma_trace = np.concatenate((gamma_trace[-(iters - burn_in_iter-1000):,:],gamma_trace_),axis=0)
					beta_trace = np.concatenate((beta_trace[-(iters - burn_in_iter-1000):,:],beta_trace_),axis=0)
					alpha_trace = np.concatenate((alpha_trace[-(iters - burn_in_iter-1000):,:],alpha_trace_),axis=0)
					top5_beta_trace = np.concatenate((top5_beta_trace[-(iters - burn_in_iter-1000):,:],top5_beta_trace_),axis = 0)

					burn_in_iter += 1000
					iters += 1000

			if (it - burn_in_iter) >= 0 and (it - burn_in_iter ) % 1000 == 0:
				print("Chain %i has sampled %i iterations " %(num,it), str(after - before),trace[it-burn_in_iter,:])

			it += 1
	
	# trace values
	trace_container[num] = {'avg': np.mean(trace,axis=0),
							'sd' : np.std(trace,axis=0)}
	del trace
	gc.collect()

	#print("trace_container[num]",num,trace_container[num]['avg'])
	#alpha values
	alpha_container[num] = {'avg': np.mean(alpha_trace,axis=0),
							'sd': np.std(alpha_trace,axis=0)}
	del alpha_trace
	gc.collect()
	#print("alpha_container[num]",num,alpha_container[num]['avg'])
	#beta values

	beta_container[num] = {'avg':np.mean(beta_trace,axis=0),
							'sd':np.std(beta_trace,axis=0)}
	del beta_trace
	gc.collect()

	gamma_container[num] = {'kmer':np.mean(gamma_trace,axis=0)}
	del gamma_trace
	gc.collect()

