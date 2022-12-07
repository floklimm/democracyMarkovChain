import numpy as _np
from collections import defaultdict as _defaultdict
from copy import deepcopy as _deepcopy
from scipy.special import loggamma as _log_gamma
from scipy.special import gamma as _gamma
from scipy.special import gammainc as _gammainc
import scipy.sparse as _sparse

from transition_counts import _k_Transitions
from transition_counts import Transitions as _Transitions
from constraint import Constraint_V2 as _Constraint

from functools import lru_cache as _lru_cache

#################################################################################################
#########################################               #########################################
#########################################     B-MON     #########################################
#########################################               #########################################
#################################################################################################

def log(x):
	return _np.log(x)
	
class BayesMultiOrderNetworkModel(object):
	"""
		This model describes which event comes after a certain chain of events.
		The chain  of events is called history `history`, while the next event is
		called `future`. Therefore the model estimates probabilities of type:
		
		p(future|history)

		The model observes paths of various length, and models the transitions
		in them. 

	"""
	def __init__(self, nodes, max_order, stationary = None, alpha_0 = 1, constraint = None):
		if constraint is None:
			constraint = _Constraint(nodes = nodes)
		self._constraint = constraint
		

		# defines the max order that we care about
		self._max_order = max_order

		# defines at which transition we assume that the paths become stationary
		if stationary is None:
			self._stationary = float("inf")
		else:
			self._stationary = stationary 
		
		self._alpha_0 = alpha_0
		self._nodes = nodes
		self._evidence = {k:_k_Transitions(k=k, 
											 nodes = nodes, 
											 constraint = constraint
											 )
							for k in range(max_order+1)}

	def observe_paths(self, paths):
		"""
			creates transitions from paths and injests them.
		"""
		transitions = _Transitions(nodes= self._nodes,
								   max_order = self._max_order,
								   constraint = self._constraint)
		transitions.add_paths(paths)
		self.observe_transitions(transitions)

	def observe_path(self, path, number_of_instances = 1):
		"""
			creates transitions using this one path and injests them.
		"""
		transitions = _Transitions(nodes = self._nodes,
								   max_order = self._max_order,
								   stationary = self._stationary,
								   constraint = self._constraint)
		transitions.add_path(path, number_of_instances = number_of_instances)
		self.observe_transitions(transitions)
	
	def observe_transitions(self, transitions):
		"""
			injests transitions by injesting each k_transition in it.
		"""
		for k in transitions.orders():
			if k<=self._max_order:
				self._observe_k_transitions(transitions.transitions_of_order(k))
			else:
				print("WARNING: Transition ignored because transition order is higher than max_order.")

	def _observe_k_transitions(self,k_transitions):
		"""
			injests _k_transitions.
		"""
		assert isinstance(k_transitions, _k_Transitions)
		self._evidence[k_transitions.get_k()] = self._evidence[k_transitions.get_k()] + k_transitions

	def loglh_paths(self, paths):
		"""
			returns the expected loglikelihood of paths
		"""
		transitions = _Transitions(nodes= self._nodes,  max_order = self._max_order, 
								  constraint = self._constraint)
		transitions.add_paths(paths)
		return self.loglh_transitions(transitions)

	def likelihood_paths(self, paths):
		return _np.exp(self.loglh_paths(paths))

	def loglh_path(self, path):
		"""
			returns the expected loglikelihood of a path
		"""
		transitions = _Transitions(nodes= self._nodes,  max_order = self._max_order, 
								  constraint = self._constraint)
		transitions.add_path(path)
		return self.loglh_transitions(transitions)

	def likelihood_path(self, path):
		return _np.exp(self.loglh_path(path))

	def loglh_transitions(self, transitions):
		""" returns the expected loglikelihood of transitions"""
		k_transitions = [transitions.transitions_of_order(order) for order in transitions.iter_orders()]
		log_L = map(self._loglh_k_transitions, k_transitions)
		return float(sum(log_L))

	def std_of_likelihood_transitions(self, transitions):
		""" returns the standard deviation of the likelihood of transitions"""
		double_transitions = _deepcopy(transitions).rescale_counters(2)
		E_sq_p = _np.exp(self.loglh_transitions(double_transitions))
		sq_E_p = _np.exp(self.loglh_transitions(transitions))**2
		return _np.sqrt( E_sq_p - sq_E_p)

	def std_of_likelihood_single_transition(self, hist, future):
		E_sq_p = _np.exp(self._loglh_per_hist(hist,{future:2}))
		sq_E_p = _np.exp(self._loglh_per_hist(hist,{future:1}))**2
		return _np.sqrt( E_sq_p - sq_E_p)

	def _loglh_k_transitions(self, k_transitions):
		"""
			This function does not do any safety checks, like:
			cheking whether everything from history is among nodes, 
			whether each possible future is in the nodes, etc.
			Do not use on its own, it's risky.
		"""
		assert k_transitions.get_k()<=self._max_order, "Order of transitions is larger than the maximum order."
		
		histories = k_transitions.get_histories()
		observed_transitions_given_history = [k_transitions.get_transitions_given_hist(hist) for hist in histories]

		log_L = map(self._loglh_per_hist, histories, observed_transitions_given_history)
		return float(sum(log_L))

	def _loglh_per_hist(self, history, observed_transitions_given_history):
		"""
			history: touple
			observed_transitions_given_history: dictionary of observed 
			transitions that we want to comute the likelihood for.
			

			The likelihood to draw x_j elements j=1..n from a multinomial
			distribution in the independent draws is the simple product:
			
				L = PRODUCT_{j = 1..n} p_j^x_j 

			Since the parameters are drawn from the Dirichlet distribution
			with hyper-parameter vector a, the expected likelihood is:
			
				E[L] = B(a+x) / B(a) 
			
			Where B(a) is multivariate Beta function

						Product{j=1..n} Gamma(a_j)
				B(a) = ----------------------------
						Gamma( Sum_{j=1..n} a_j )
			
			Btw, this is a fraction, xD.

			I compute the log of this function, which means the product 
			becomes a sum.
			
				E[logL]	= 	+ sum(loggamma(a_j+x_j)) 
							- loggamma(sum(a_j+x_j)) 
					  		- sum(loggamma(a_j))
					  		+ loggamma(sum(a_j)) 
		"""
		remaining_futures = self._constraint.future_of(history)

		for future in observed_transitions_given_history:
			if future not in remaining_futures:
				return -float("inf")

		sum_a = 0
		sum_loggamma_ai = 0
		sum_x = 0
		sum_loggamma_ai_xi = 0
		
		for future in observed_transitions_given_history:
			a_i = self._get_alpha(history, future)
			
			sum_a += a_i
			sum_loggamma_ai += _log_gamma(a_i)

			x_i = observed_transitions_given_history[future]
			sum_x += x_i
			sum_loggamma_ai_xi += _log_gamma(a_i + x_i)
			remaining_futures.remove(future)

		for future in remaining_futures:
			a_i = self._get_alpha(history, future)
			sum_a += a_i

		return float(sum_loggamma_ai_xi - _log_gamma(sum_a+sum_x) - sum_loggamma_ai + _log_gamma(sum_a))

	def _get_alpha(self, history, future):
		assert len(history) <= self._max_order, "History {} too long. It can be at most {} long.".format(history, self._max_order)
		n_path = self._evidence[len(history)].get_transition_counter(history, future)
		alpha = self._alpha_0 + n_path
		assert alpha >= 0, "Alpha should be non-negative"
		return alpha

	def adjacency_matrix(self, order):
		assert order<=self._max_order, "Order has to be at most max_order. order = "+str(order)
		assert order>0, "Order has to be larger than 0. order = "+str(order)		

		hon_node_to_index = defaultdict(lambda:-1)
		index_to_hon_node = {}
		index = 0
		def update_maps(hist):
			if hon_node_to_index[hist] == -1:
				hon_node_to_index[hist] = index
				index_to_hon_node[index] = hist
				index += 1

		node_to = []
		node_from = []
		value = []
		for hist in self._constraint.iterate_possible_histories(order):
			for future in self._constraint.future_of(hist):
				node_2 = (*hist[1:], future)
				update_maps(hist)
				update_maps(node_2)
				node_from.append(hon_node_to_index[hist])
				node_to.append(hon_node_to_index[node_2])
		
		shape = (index+1, index+1)
		value = _np.repeat(1, len(node_from))
		A = _sparse.coo_matrix((value, (node_from, node_to)), shape=shape).tocsr()
		return A, hon_node_to_index, index_to_hon_node

	def transition_matrix(self, order):
		assert order<=self._max_order, "Order has to be at most max_order. order = "+str(order)
		assert order>0, "Order has to be larger than 0. order = "+str(order)		

		hon_node_to_index = defaultdict(lambda:-1)
		index_to_hon_node = {}
		index = 0
		def update_maps(hist):
			if hon_node_to_index[hist] == -1:
				hon_node_to_index[hist] = index
				index_to_hon_node[index] = hist
				index += 1

		node_to = []
		node_from = []
		value = []
		for hist in self._constraint.iterate_possible_histories(order):
			for future in self._constraint.future_of(hist):
				node_2 = (*hist[1:], future)
				update_maps(hist)
				update_maps(node_2)
				node_from.append(hon_node_to_index[hist])
				node_to.append(hon_node_to_index[node_2])
				value.append( self._loglh_per_hist(hist, {future:1}))
		shape = (index+1, index+1)
		T = _sparse.coo_matrix((value, (node_from, node_to)), shape=shape).tocsr()
		return T, hon_node_to_index, index_to_hon_node

	def laplacian_matrix(self,):
		L = ...
		return L

	def plot_explanation_path_loglikelihood(self, path, figsize = (10,10), filter_histories = None):
		transitions = _Transitions(nodes= self._nodes,  max_order = self._max_order, 
								  constraint = self._constraint)
		transitions.add_path(path)
		self.plot_explanation_transitions_loglikelihood(transitions, figsize = figsize, filter_histories = filter_histories)
	
	def plot_explanation_paths_loglikelihood(self, paths, figsize = (10,10), filter_histories = None):
		transitions = _Transitions(nodes= self._nodes,  max_order = self._max_order, 
								  constraint = self._constraint)
		transitions.add_paths(paths)
		self.plot_explanation_transitions_loglikelihood(transitions, figsize = figsize, filter_histories = filter_histories)

	def plot_explanation_transitions_loglikelihood(self, transitions, figsize = (10,10), filter_histories = None):
		import matplotlib.pyplot as plt

		for k in transitions.iter_orders():
			k_transitions = transitions.transitions_of_order(k) 

			histories = k_transitions.get_histories()
			if filter_histories is None:
				pass
			elif isinstance(filter_histories, int):
				histories = [hist for hist in histories if k_transitions.get_history_counter(hist)>filter_histories]
			else:
				histories = filter_histories(histories)

			histories.sort(key = lambda hist: -k_transitions.get_history_counter(hist))
			
			n_hist = len(histories)
			
			# scale the numbers of plots per x and y axis
			# so that the plots look more or less square.
			if figsize[1]<figsize[0]:
				product_times_ratio = n_hist *figsize[1]/figsize[0]
				n = int(_np.ceil(_np.sqrt(product_times_ratio))) 
				m = int(_np.ceil(int(n_hist/n))) 
			else:
				product_times_ratio = n_hist *figsize[0]/figsize[1]
				m = int(_np.ceil(_np.sqrt(product_times_ratio)) )
				n = int(_np.ceil(n_hist/m) )

		
			ij_from_hist = {hist:(int(e/m), e%m) for e,hist in enumerate(histories)}
			fig, ax = plt.subplots(n, m, figsize = figsize)
			
			# the ax should always be 2D array
			if n == 1:
				ax = [ax]
			if m == 1:
				ax = [[axi] for axi in ax]

			for hist in histories:
				i,j = ij_from_hist[hist]
				futures = [xi for xi in self._constraint.future_of(hist)]
				futures.sort(key = lambda x: - k_transitions.get_transition_counter(hist, x))
				x = [str(xi) for xi in futures]
				obs_y = _np.array([k_transitions.get_transition_counter(hist, xi) for xi in x])
				obs_y = obs_y/_np.sum(obs_y)
				exp_y = _np.array([self._get_alpha(hist, xi) for xi in x])
				exp_y = exp_y/_np.sum(exp_y)
				
				std_y = [self.std_of_likelihood_single_transition(hist, future) for future in futures]

				def map_between_01(x):
					if x>40:
						return 1
					else:
						return _np.exp(x)/(100+_np.exp(x))
				alpha_overall = 0.2+0.8*(map_between_01(k_transitions.get_history_counter(hist)))
				ax[i][j].bar(x, height = obs_y, color = (0,0.5,0), alpha = 0.5*alpha_overall)
				ax[i][j].bar(x, height = exp_y, yerr = std_y,ecolor  = (0.7,0.35,0.35), color = (0.5,0,0), alpha = 0.5*alpha_overall)
				ax[i][j].set_xticks(x, x)
				ax[i][j].set_title(hist)
				ax[i][j].legend([],title = "{} observations".format(k_transitions.get_history_counter(hist)) )
			plt.show()
	
	def generate_random_path_of_length(self, l):
		path = tuple()
		
		@_lru_cache(maxsize = None)
		def futures_and_probabilities(hist):
			futures = tuple(self._constraint.future_of(hist))
			alphas = tuple(self._get_alpha(hist, future) for future in futures)
			probabilities = tuple(_np.random.dirichlet(alphas))
			return (futures, probabilities)

		for i in range(l+1):
			hist = path[-self._max_order:]
			futures, probs = futures_and_probabilities(hist)
			# print(futures)
			# print(probs)
			# print(hist)
			next_node = _np.random.choice(futures, p =  probs)
			path = path + (next_node,)
		futures_and_probabilities.cache_clear()
		return path

	def _differential_entropy_of_hist(self,history):
		"""
			history: touple
			
			returns: differential entropy of the parameters p,
					 that we use to model p(future_i|hist)

			Note: this entropy can be negative. And is weird, generally.

			Formula for entropy is:
				H[P] = log B(a) + (a - n) digamma(a) - SUM_{j = 1...n} (a_j - 1)*digamma(aj)

			Where B(a) is multivariate Beta function

						Product{j=1..n} Gamma(a_j)
				B(a) = ----------------------------
						Gamma( Sum_{j=1..n} a_j )
			
			Btw, this is a fraction, xD.

			I compute the log of this function, which means the product 
			becomes a sum.
			
				log B(a) =  sum( loggamma(a_j)) - loggamma(sum(a_j)) 
			
		"""
		futures = self._constraint.future_of(history)
		n = len(futures)

		log_B_a = 0
		sum_a = 0
		sum_ai_digamma_ai = 0
		for future in futures:
			a_i = self._get_alpha(history, future)
			
			sum_a += a_i
			log_B_a += _log_gamma(a_i)
			
			sum_ai_digamma_ai += (a_i - 1)*_digamma(a_i)
		
		log_B_a -= _log_gamma(sum_a)

		H = log_B_a + (sum_a - n)*_digamma(sum_a) - sum_ai_digamma_ai
		return H

	# def model_guess_for_path_end(self, path_beginning, threshold):
		# """
		# 	path_beginning is the known part of the path.
		# 	threshold is the probability.

		# 	method returns all endings of the path that have probability higher than this threshold.
		# """
		# history = path_beginning
		# for future in self._constraint.future_of(history):
		# 	if 

def ML_loglh_and_DoF_of_BMONs(paths, max_order, nodes, stationary = None, constraint = None):
	if constraint is None:
		constraint = _Constraint(nodes)

	log_lh = {}
	dof = {}
	opt_order = 0
	
	n_obs = _defaultdict(int)
	p_vals = []
	for order in range(0,max_order+1):
		transitions = _Transitions(nodes,  order, stationary = stationary, constraint = constraint)

		transitions.add_paths(paths)

		llh = 0
		
		for leng in paths.lengths():
			for path in paths.paths_of_length(leng):
				for e, node in enumerate(path):
					n_obs[order] += paths.path_counter(path)

					future = node
					hist = path[max(0,e-order):e]
					hist_counter = transitions.get_history_counter(hist)
					future_counter = transitions.get_transitions_given_hist(hist)[future]
					llh += log(future_counter/hist_counter) * paths.path_counter(path) 
					

		log_lh[order] = llh

		dof[order] = constraint.degrees_of_freedom(order)
	return log_lh, dof, n_obs

def AIC(ln_lh, dof):
	return 2*dof - 2*ln_lh


def BIC(ln_lh, dof, n_obs):
	return _np.log(n_obs)*dof - 2*ln_lh


def EDC(ln_lh, dof, n_obs, n_nodes):
	return _np.log(_np.log(n_obs))*dof/(n_nodes-1) - 2*ln_lh
	
def opt_order_likelihood_ratio_test(p_level, paths, max_order, nodes, stationary = None, constraint = None, alpha_0 = 1, noprint = True):
	if constraint is None:
		constraint = _Constraint(nodes)

	log_lh = {}
	dof = {}
	opt_order = 0
	
	p_vals = []
	for order in range(0,max_order+1):
		transitions = _Transitions(nodes,  order, stationary = stationary, constraint = constraint)

		transitions.add_paths(paths)

		llh = 0
		
		for leng in paths.lengths():
			for path in paths.paths_of_length(leng):
				for e, node in enumerate(path):
					future = node
					hist = path[max(0,e-order):e]
					hist_counter = transitions.get_history_counter(hist)
					future_counter = transitions.get_transitions_given_hist(hist)[future]
					llh += log(future_counter/hist_counter) * paths.path_counter(path) 

		log_lh[order] = llh

		dof[order] = constraint.degrees_of_freedom(order)

		
		if order>0:
			p_val = p_value(log_lh[order-1], log_lh[order], dof[order-1], dof[order])
			p_vals.append(p_val)
			reject = (p_level>p_val)
			if reject:
				if opt_order == order -1:
					opt_order = order
	
	if not noprint:
		print("log_lh:",log_lh)
		print("dof:   ",dof)
	return opt_order, p_vals

def p_value(loglikelihood1, loglikelihood2, dof1, dof2):
	degree_diff = dof2-dof1
	minus_two_log_ratio_lh =  - 2 * (loglikelihood1 - loglikelihood2)
	p = 1 - cummulative_chi_sqare(degree_diff, minus_two_log_ratio_lh)
	return p

def cummulative_chi_sqare(k, x):
	# it should be incomplete_gamma(k,x)/gamma(k)
	# but gammainc is already normalied with gamma(k)
	return _gammainc(k/2, x/2)

def opt_order_bayes_factors(paths, max_order, nodes, eps = 1, stationary = None, constraint = None, alpha_0 = 1):
	"""
		opt_order is the largest order that is eps times more likely 
		than every order smaller than it. Default eps = 1 means the 
		one that is the most likely is the optimal one.
	"""

	log_lh = {}
	for order in range(0,max_order+1):
		model = BayesMultiOrderNetworkModel(nodes, order,stationary = stationary, alpha_0 = alpha_0, constraint = constraint)
		log_lh[order] = model.loglh_paths(paths)
	
	opt_order = 0
	for order in range(1, max_order+1):
		opt = True
		for alt_order in range(order):
			if log_lh[order] < log_lh[alt_order] + log(eps):
				opt = False
		if opt:
			opt_order = order

	return opt_order, log_lh

def probability_of_orders(log_lh):
	"""
		returns the dictionary of probabilities for each order.
	"""
	normalize = max(log_lh.values())
	def loglikelihood_to_prob(loglikelihood):
		if loglikelihood - normalize < -700:
			return 0
		else:
			return _np.exp(loglikelihood - normalize)
	p_of_k = {k:loglikelihood_to_prob(log_lh[k]) for k in log_lh}
	suma = sum(p_of_k.values())
	for k in p_of_k:
		p_of_k[k] /= suma
	return p_of_k