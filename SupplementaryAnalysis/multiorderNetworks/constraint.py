from copy import deepcopy
from collections import defaultdict as _defaultdict

class Constraint():
	"""
		The function of this class is to define when a future could follow a history.
	"""
	def __init__(self, nodes, path_end = "end", model_path_ends = False):
		self._nodes = set(nodes)
		assert len(self._nodes) > 0, "Nodes cannot be empty."
		
		self._model_path_ends = model_path_ends
		if model_path_ends:
			self._path_end = path_end

		self._constraints = []

	def check(self, hist, future):
		"""
			Apply _constraints on the history/future pair.
		"""
		assert self._basic_transition(hist, future)
		for constraint in self._constraints:
			if not constraint(hist, future):
				return False
		return True

	def future_of(self, hist):
		"""
			returns list of possible futures.
		"""
		possibilities = deepcopy(self._nodes)
		if self._model_path_ends:
			possibilities.add(self._path_end)
		for c in self._constraints:
			possibilities = [x for x in possibilities if c(hist, x)]
		return possibilities

	def degrees_of_freedom(self, order):
		"""
			computes the degrees of freedom 
			for a model of order up to `order`
		"""
		return self._degrees_of_freedom_given_hist(tuple(), order)

	def _degrees_of_freedom_given_hist(self, hist, order):
		"""
			computes (recursively) the degrees 
			of freedom for a model of order up 
			to `order` that govern the future of
			a history `hist`.
		"""
		possibilities = self.future_of(hist)
		if possibilities != []:
			dof = len(possibilities) - 1
			if order == len(hist):
				return dof 
			else:
				return dof + sum([self._degrees_of_freedom_given_hist((*hist, node), order) 
								  for node in possibilities])
		else:
			return 0

	def _basic_transition(self, hist, future):
		"""
			returns error if:
				1. _path_end is in the history.
				2. future is neither _path_end nor in the _nodes.
		"""
		assert isinstance(hist, tuple), "Parameter hist has to be an instance of tuple."
		assert all([x in self._nodes for x in hist]), "Encountered node in history which is not an element of nodes"
		if self._model_path_ends:
			assert future in self._nodes or future == self._path_end, "Future {} is not an element of nodes and not _path_end".format(future, self._nodes)
		else:
			assert future in self._nodes, "Future {} is not an element of nodes.".format(future)
		return True

	def delete_added_constraints(self):
		"""
			resets the _constraints
		"""
		self._constraints = []

	def add_network_constraint(self, edges):
		"""
			network constraint.
		"""

		assert all([x in self._nodes and y in self._nodes for x,y in edges]), "Unknown nodes."

		if self._model_path_ends:
			def network_constraint(hist, future):
				if future == self._path_end:
					return True
				elif hist ==tuple():
					return (future in self._nodes)
				else:
					# if it is not the beginning or the end, the link from history to the future has to exist
					return (future in self._nodes) and ((hist[-1],future) in edges)
		else:
			def network_constraint(hist, future):
				future_ok = future in self._nodes
				if hist != tuple():
					hist_ok = ((hist[-1],future) in edges)
				else:
					hist_ok = True

				return future_ok and hist_ok
		# set this function as a constraint
		self._constraints.append(network_constraint)

	def add_no_backtracking(self):
		"""
			Introduces no-backtracking constraint.
			This means that the future cannot appear in the history.
			Note that if the history is not the full hisotry, the 
			no-backtracking constraint will be checked only for the given history
		"""
		def no_backtracking_constraint(hist, future):
			return future not in hist
		self._constraints.append(no_backtracking_constraint)
	
	def add_custom_constraint(self,f):
		"""
			introduces constraint on the underlying topology.
			f is the constraint:
				input is hist/future pair, 
				returns boolean value.
			If f(hist, future) = False then the prior of p(future|hist) = 0.
			if it is not introduced, the constraint allow everything
		"""
		assert isinstance(f((list(self._nodes)[0],), list(self._nodes)[0]), bool), "f should take history/future pair and return a boolean value"
		self._constraints.append(f)

	def iterate_possible_histories(self,order):
		assert False, "Not implemented yet."

	def prolong_the_paths_by_one_step(self, paths):
		"""
			Output is the list of paths that is obtained by prolonging every path from the set `paths`.
		"""
		new_paths = set()
		for path in paths:
			for future in self.future_of(path):
				new_paths.add(path+(future,))
		return new_paths

	def all_paths_of_length(self, length, up_to=False):
		"""
			returns all possible paths of order (up to) some length
		"""
		assert length > 0

		if not up_to:
			paths = self.prolong_one_step([tuple()])
			#now they are of length 0. (only starting nodes)
			for i in range(length):
				paths = prolong_one_step(paths)
			# now they are of proper length
		else:
			l_p_c = self.prolong_one_step(paths)
			paths = self.prolong_one_step(paths)
			# now they are of length 0.(only starting nodes)
			for i in range(length):
				l_p_c = self.prolong_one_step(l_p_c)
				paths = paths + l_p_c
			# now they are of proper length

		return paths
		
	# def relabel_nodes(self, node_to_label):
	# 	assert False, "Not implemented yet."

	def __eq__():
		"""
			This function is needed to compare two _constraints.
			For example, when I add two transitions objects.
		"""
		assert False, "Not implemented yet."


class Constraint_V2():
	"""
		The function of this class is to define when a future could follow a history.
	"""
	def __init__(self, nodes, mode = "safe"):
		"""
			mode can be "performance" or "safe". Safe mode will check every time 
			whether all nodes mentioned in history or future are among nodes.
			Performance mode will skip this test.

		"""
		assert len(nodes) > 0, "Nodes cannot be empty."
		self._nodes = set(nodes)
		
		assert mode in ["safe", "performance"]
		self._mode = mode

	def check(self, hist, future):
		"""
			Check the history/future pair.
		"""
		assert self._basic_transition(hist, future)
		return True

	def future_of(self, hist):
		"""
			returns list of possible futures.
		"""
		assert self._valid_history(hist)
		# We use the fact that it is a list in degrees of freedom.
		return list(self._nodes)
		
	def degrees_of_freedom(self, order):
		"""
			computes the degrees of freedom 
			for a MON model of max_order equal `order`
		"""
		assert isinstance(order, int) and order >= 0, "Order has to be non-negative integer."	
		if type(self) == Constraint_V2:
			N = len(self._nodes)
			return N**order * (N-1)
		else:
			return self._degrees_of_freedom_given_hist(tuple(), order)

	def _degrees_of_freedom_given_hist(self, hist, order):
		"""
			computes (recursively) the degrees 
			of freedom for a model of order up 
			to `order` that govern the future of
			a history `hist`.
		"""
		possibilities = self.future_of(hist)
		if possibilities != []:
			dof = len(possibilities) - 1
			if order == len(hist):
				return dof 
			else:
				return dof + sum([self._degrees_of_freedom_given_hist((*hist, node), order) 
								  for node in possibilities])
		else:
			return 0

	def _valid_history(self, hist):
		if self._mode == "safe":
			assert isinstance(hist, tuple), "Parameter hist has to be an instance of tuple."
			assert all([x in self._nodes for x in hist]), "Encountered node in history {} which is not an element of nodes \n {}".format(hist, self._nodes ) 
		return True

	def _valid_future(self, future):
		if self._mode == "safe":
			assert future in self._nodes, "Future {} is not an element of nodes.".format(future)
		return True

	def _basic_transition(self, hist, future):
		"""
			In safe mode, returns error if:
				1. _path_end is in the history.
				2. future is neither _path_end nor in the _nodes.
		"""
		assert self._valid_future(future)
		assert self._valid_history(hist)
		return True

	def prolong_the_paths_by_one_step(self, paths):
		"""
			Output is the list of paths that is obtained by prolonging every path from the set `paths`.
		"""
		new_paths = list()
		for path in paths:
			for future in self.future_of(path):
				new_paths.append(path+(future,))
		return new_paths

	def all_paths_of_length(self, length, up_to=False):
		"""
			returns all possible paths of order (up to) some length
		"""
		assert length > 0

		if not up_to:
			paths = self.prolong_the_paths_by_one_step([tuple()])
			#now they are of length 0. (only starting nodes)
			for i in range(length):
				paths = self.prolong_the_paths_by_one_step(paths)
			# now they are of proper length
		else:
			l_p_c = self.prolong_the_paths_by_one_step(paths)
			paths = self.prolong_the_paths_by_one_step(paths)
			# now they are of length 0.(only starting nodes)
			for i in range(length):
				l_p_c = self.prolong_the_paths_by_one_step(l_p_c)
				paths = paths + l_p_c
			# now they are of proper length

		return paths

class BagOfConstraints(Constraint_V2):
	"""
		The function of this class is to define when a future could follow a history.
	"""
	def __init__(self, nodes, mode = "safe"):
		"""
			mode can be "performance" or "safe". Safe mode will check every time 
			whether all nodes mentioned in history or future are among nodes.
			Performance mode will skip this test.

		"""
		Constraint_V2.__init__(self, nodes, mode = mode)

		self._constraints = []

	def check(self, hist, future):
		"""
			Apply _constraints on the history/future pair.
		"""
		assert self._basic_transition(hist, future)
		
		for constraint in self._constraints:
			if not constraint.check(hist, future):
				return False
		return True

	def future_of(self, hist):
		"""
			returns list of possible futures.
		"""
		assert self._valid_history(hist)
		if len(self._constraints) == 0:
			return self._nodes
		else:
			# compute the intersection of all
			possibilities = self._constraints[0].future_of(hist)
			for c in self._constraints[1:]:
				new_possibilities = []
				allowed = c.future_of(hist)
				# keep only those that appear in both.
				for i in possibilities:
					if i in allowed:
						new_possibilities.append(i)
				possibilities = new_possibilities
			return possibilities

	def delete_constraints(self):
		"""
			resets the _constraints
		"""
		self._constraints = []

	def add_constraint(self, constraint):
		assert isinstance(constraint, Constraint_V2) 
		self._constraints.append(constraint)

class NetworkConstraint(Constraint_V2):
	"""Network Constraint Class"""
	def __init__(self, nodes, edges, mode = "safe"):
		Constraint_V2.__init__(self, nodes, mode = mode)

		self._successors = _defaultdict(set)
		for s, t in edges:
			assert self._basic_transition((s,), t)
			self._successors[s].add(t)

	def check(self, hist, future):
		assert self._basic_transition(hist, future)
		if hist == tuple():
			return True
		else:
			return future in self._successors[hist[-1]]

	def future_of(self, hist):
		assert self._valid_history(hist)
		if hist == tuple():
			return list(self._nodes)
		else:
			return list(self._successors[hist[-1]])
		
class NoBacktrackingConstraint(Constraint_V2):
	"""Forbids Backtracking"""
	def __init__(self, nodes, mode = "safe"):
		Constraint_V2.__init__(self, nodes, mode)

	def check(self, hist, future):
		"""
			Check the history/future pair.
		"""
		assert self._basic_transition(hist, future)
		return future not in hist
	
	def future_of(self, hist):
		"""
			returns list of possible futures.
		"""
		assert self._valid_history(hist)	
		return list(self._nodes - set(hist))
		