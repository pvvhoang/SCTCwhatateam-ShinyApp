#!/usr/bin/python -tt

#' @export

import numpy as np
from itertools import cycle
import pandas
# Use unsupecised stepwise feature selection implemented in https://github.com/EFavDB/linselect

class Base(object):
	"""
	Base class for the stepwise selection algorithms.

	Attributes common to all child classes
	--------------------------------------
	M : np.array (self.dimension, self.dimension)
		The correlation matrix of all variables.  Numpy automatically scales
		all variables so that each has mean 0 and variance 1, as needed.

	C : np.array (self.dimension, self.dimension)
		This holds all information needed to identify the gains and costs
		associated with moving features into or out of the predictor set at
		each step of the process.

	s : np.array, Boolean (self.dimension, )
		This array's index i specifies whether variable i is currently in the
		predictor set.

	mobile : np.array, Boolean (self.dimension, )
		This array's index i specifies whether variable i is free to move in
		and out of the predictor set.

	targets : np.array, Boolean (self.dimension, )
		This array's index i specifies whether variable i is in the target set,
		i.e., one of the variables we are trying to fit.

	dimension : int
		This is the number of variables, including both the features in X and
		the targets in y (when passed).

	dtype : numeric variable type
		Computations will be carried out using this level of precision.
		Note: Lower precision types result in faster computation. However, for
		nearly redundant data sets these can sometimes result in nan results
		populating the ordered_cods.

	Private methods
	---------------
	_setup(self, X, s, mobile, targets)
		Sets passed search parameters.

	_set_efroymson(self)
		Sets the Efroymson matrix C.

	_set_cod(self)
		Evaluates current COD of the target set.

	_forward_step(self)
		Identifies optimal forward step and COD gain associated with this.

	_reverse_step(self)
		Identifies optimal reverse step and COD cost associated with this.
	"""

	def __init__(self, dtype):
		self.dtype = dtype
		self.M = None
		self.s = None
		self.mobile = None
		self.targets = None

	def _setup(self, X=None, s=None, mobile=None, targets=None):
		"""
		Reset various parameters if passed.  If not passed and values are not
		currently set, default to forward selection initial placements.
		"""
		if X is None and self.M is None:
			raise ValueError('Must pass X to initialize')
		if X is not None:
			if self.M is not None:
				raise ValueError('X was previously passed')
			X = np.array(X)
			self.dimension = X.shape[1]
			self.indices = np.arange(self.dimension)
			self.M = np.corrcoef(np.transpose(X)).astype(self.dtype)
		if s is not None:
			self.s = np.array(s)
		elif self.s is None:
			self.s = np.full(self.dimension, False, dtype=bool)
		if mobile is not None:
			self.mobile = np.array(mobile)
		elif self.mobile is None:
			self.mobile = np.full(self.dimension, True, dtype=bool)
		if targets is not None:
			if self.targets is not None:
				raise ValueError('targets was previously set')
			self.targets = np.array(targets)
		elif self.targets is None:
			self.targets = np.full(self.dimension, True, dtype=bool)

	def _set_efroymson(self):
		"""
		Set the Efroymson matrix, C.
		"""
		self.C = (np.diag(~self.s).dot(self.M) + np.diag(self.s)).dot(
			np.linalg.inv(np.diag(self.s).dot(self.M).dot(np.diag(self.s))
						  + np.diag(~self.s)) - np.diag(~self.s)).dot(
			self.M.dot(np.diag(~self.s)) + np.diag(self.s))

	def _set_cod(self):
		"""
		Evaluate current COD of target set.
		"""
		self.cod = (np.diag(self.C)[~self.s * self.targets].sum()
			+ self.dtype(np.sum(self.s * self.targets)))

	def _forward_step(self):
		"""
		Take optimal forward step and update C, s, and COD.  Return optimal
		index.
		"""
		# check there are mobile options, exit if not
		if np.sum(self.mobile * ~self.s) == 0:
			return

		# identify COD gain of each candidate move, select best
		gains = np.einsum('ij,ij->i',
			(self.C - self.M)[np.ix_(
				~self.s * self.mobile, ~self.s * self.targets)],
			(self.C - self.M)[np.ix_(
				~self.s * self.mobile, ~self.s * self.targets)]) / (
				self.dtype(1.0) - self.C.diagonal()[np.ix_(
					~self.s * self.mobile)])
		opt_index = np.argmax(gains)
		opt_gain = gains[opt_index]

		# map opt_index to index in full feature set
		opt_index = self.indices[np.ix_(
			~self.s * self.mobile)][opt_index]

		# save what the moved column should change to
		C_opt_index_update = (
			self.M[opt_index] * ~self.s - self.C[opt_index]) / (
				self.dtype(1.0) - self.C[opt_index, opt_index])
		C_opt_index_update[opt_index] = self.dtype(1.0) / (
			self.dtype(1.0) - self.C[opt_index, opt_index])

		# update C
		x = self.C[opt_index] - self.M[opt_index] * ~self.s
		for index in range(self.dimension):
			self.C[index] -= x[index] * x / x[opt_index]

		# fix the opt_index column and row
		self.C[opt_index, :] = C_opt_index_update
		self.C[:, opt_index] = C_opt_index_update

		# update s
		self.s[opt_index] = True

		# update COD, refresh if update gives nan
		test_cod = self.cod + opt_gain
		if np.isnan(test_cod):
			self._set_cod()
		else:
			self.cod = test_cod

		# return opt_index
		return opt_index

	def _reverse_step(self):
		"""
		Take optimal reverse step and update C, s, and COD.  Return optimal
		index.
		"""
		# check there are mobile options, exit otherwise
		if np.sum(self.mobile * self.s) == 0:
			return

		# evaluate COD costs
		costs = (
			self.targets[self.s * self.mobile]
			+ np.einsum('ij,ij->i',
				self.C[np.ix_(
					self.s * self.mobile, ~self.s * self.targets)],
				self.C[np.ix_(
					self.s * self.mobile, ~self.s * self.targets)])
			) / (self.C.diagonal()[self.s * self.mobile])
		opt_index = np.argmin(costs)
		opt_cost = costs[opt_index]

		# map opt_index to index in full feature set
		opt_index = self.indices[np.ix_(self.s * self.mobile)][opt_index]

		# save what the moved column should change to
		C_opt_index_update = (
			self.M[opt_index] * ~self.s - self.C[opt_index] / (
			self.C[opt_index, opt_index]))
		C_opt_index_update[opt_index] = self.M[opt_index, opt_index] - (
			self.dtype(1.0) / self.C[opt_index, opt_index])

		# update C
		x = self.C[opt_index].copy()
		for index in range(self.dimension):
			self.C[index] -= x[index] * x / x[opt_index]

		# fix the top_index column and row
		self.C[opt_index, :] = C_opt_index_update
		self.C[:, opt_index] = C_opt_index_update

		# update s
		self.s[opt_index] = False

		# update COD, refresh if update gives nan
		test_cod = self.cod - opt_cost
		if np.isnan(test_cod):
			self._set_cod()
		else:
			self.cod = test_cod

		# return opt_index
		return opt_index
	
class FwdSelect(Base):
	"""
	FwdSelect -- Efficient Forward Stepwise Linear Regression

	A class for carrying out forward, single-step linear feature selection
	protocols.  At each step, the feature that is selected is that which
	increases the total COD (coefficient of determination, aka R^2) by the
	largest amount.  The feature ordering and CODs are stored, allowing for
	review.

	Special Attributes
	------------------
	ordered_features: list
		List of the feature indices.  The ordering is that in which the
		features were added to the predictor set during selection.

	ordered_cods: list
		This list's index i specifies the COD that results if only the first i
		features of ordered_features are taken as predictors (large COD
		values are better and a perfect score = n_targets).
	"""
	def __init__(self, dtype=np.float32):
		"""
		input
		-----
		dtype: numeric variable type
		"""
		super(FwdSelect, self).__init__(dtype)

	def fit(self, X, y=None):
		"""
		Method fits passed data, evaluates self.ordered_features and
		self.ordered_cods.

		Parameters
		----------
		X : np.array (n_examples, n_features)
			Data array of features containing samples across all features. Must
			be numeric.

		y : np.array (n_examples, n_targets), default None
			Array of label values for each example. If n_targets > 1 we seek
			the features that maximize the sum total COD over the separate
			labels.  If None passed, we carry out unsupervised selection,
			treating all features as targets.  If passed, must be numeric.

		Returns
		-------
		self : fitted instance
		"""
		# setup
		if y is not None:
			X = np.hstack((X, y))
			self._setup(X)
			self.mobile[-y.shape[1]:] = False
			self.targets = ~self.mobile
		else:
			self._setup(X)
		self._set_efroymson()
		self._set_cod()
		self.ordered_cods = list()
		self.ordered_features = list()

		# carry out stepwise procedure
		for _ in range(np.sum(self.mobile)):
			opt_index = self._forward_step()
			self.ordered_features.append(opt_index)
			self.ordered_cods.append(self.cod)
		return self
	
class RevSelect(Base):
	"""
	RevSelect -- Efficient Reverse Stepwise Linear Regression

	A class for carrying out reverse, single-step linear feature selection
	protocols.  At each step, the feature that is selected is that which
	reduces the total COD (coefficient of determination, aka R^2) by the least
	amount.  The feature ordering and CODs are stored, allowing for review.

	Special Attributes
	------------------
	ordered_features: list
		List of the feature indices.  The ordering is the reverse of that in
		which the features were removed from the predictor set during
		selection.

	ordered_cods: list
		This list's index i specifies the COD that results if only the first i
		features of ordered_features are taken as predictors (large COD
		values are better and a perfect score = n_targets).
	"""
	def __init__(self, dtype=np.float32):
		"""
		Parameters
		----------
		dtype: numeric variable type
			Computations will be carried out using this level of precision.
			Note: Lower precision types result in faster computation. However,
			for nearly redundant data sets these can sometimes result in nan
			results populating the ordered_cods.
		"""
		super(RevSelect, self).__init__(dtype)

	def fit(self, X, y=None):
		"""
		Method fits passed data, evaluates self.ordered_features and
		self.ordered_cods.

		Parameters
		----------
		X : np.array (n_examples, n_features)
			Data array of features containing samples across all features. Must
			be numeric.

		y : np.array (n_examples, n_targets), default None
			Array of label values for each example. If n_targets > 1 we seek
			the features that maximize the sum total COD over the separate
			labels.  If None passed, we carry out unsupervised selection,
			treating all features as targets.  If passed, must be numeric.

		Returns
		-------
		self : fitted instance
		"""
		# setup
		if y is not None:
			X = np.hstack((X, y))
			self._setup(X)
			self.mobile[-y.shape[1]:] = False
			self.targets = ~self.mobile
		else:
			self._setup(X)
		self.s = self.mobile.copy()
		self._set_efroymson()
		self._set_cod()
		self.ordered_cods = [self.cod]
		self.ordered_features = list()
		# carry out stepwise procedure
		for _ in range(np.sum(self.mobile)):
			opt_index = self._reverse_step()
			self.ordered_features.append(opt_index)
			self.ordered_cods.append(self.cod)
		# reverse list orders for consistency with FwdSelect
		self.ordered_features = self.ordered_features[::-1]
		self.ordered_cods = self.ordered_cods[-2::-1]
		return self
	
class GenSelect(Base):
	"""
	GenSelect -- Efficient General Stepwise Linear Regression

	A class for carrying out general, single-step linear feature selection
	protocols: protocols that include both forward and reverse search steps.
	Best results seen so far are stored, allowing for review.  This also allows
	for a search to be continued with repositioning or step protocol
	adjustments, as desired.

	Special Attributes
	------------------
	best_results: dict
		Keys of this dict correspond to feature subset size.  The value for a
		given key is also a dict -- one characterizing the best subset seen so
		far of this size.  These inner dicts have two keys, `s` and `cod`.  The
		first holds a Boolean array specifying which features were included in
		the subset and the second holds the corresponding COD.
	"""
	def __init__(self, dtype=np.float32):
		"""
		input
		-----
		dtype: numeric variable type
			Computations will be carried out using this level of precision.
			Note: Lower precision types can result in faster computation.
			However, for nearly redundant data sets these can sometimes result
			in nan results populating the ordered_cods.
		"""
		super(GenSelect, self).__init__(dtype)
		self.best_results = {}

	def _update_best_results(self):
		"""
		If current COD is larger than prior best at current predictor set size,
		update best_results dict.
		"""
		# do not update if current COD is nan
		if np.isnan(self.cod):
			return
		# possible update
		size_s = np.sum(self.s)
		if size_s in self.best_results.keys():
			if self.cod <= self.best_results[size_s]["cod"]:
				# better result before, do not update
				return
		new_dict = {"s": self.s.copy(), "cod": self.cod}
		self.best_results[size_s] = new_dict

	def position(self, X=None, s=None, mobile=None, targets=None):
		"""
		Set the operating conditions of the stepwise search.

		Parameters
		----------
		X: np.array, (n_examples, n_features), default None
			The data set to be fit.  This must be passed in the first call to
			this method, but should not need to be passed again in any
			following repositioning call.  Must be numeric.

		s: np.array, (n_features), default None
			This is a Boolean array that specifies which predictor set to use
			when we begin (or continue) our search.  If the index i is set to
			True, the corresponding feature i will be included in the initial
			predictor set.

		mobile: np.array, (n_features), default None
			This is a Boolean array that specifies which of the features are
			locked into or out of our fit -- if the index i is set to True, the
			corresponding feature i can be moved into or out of the predictor
			set.  Otherwise, the feature i is locked in the set specified by
			the passed s argument.

		targets: np.array, (n_features), default None
			This is a Boolean array that specifies which of the columns of X
			are to be fit -- analogs of y in the FwdSelect and RevSelect
			algorithms.  If the index i is set to True, the corresponding
			column i will be placed in the target set.  Once set, this should
			not be passed again in any following repositioning call.
		"""
		self._setup(X, s, mobile, targets)
		self._set_efroymson()
		self._set_cod()
		self._update_best_results()

	def search(self, protocol=(2, 1), steps=1):
		"""
		(Continue) stepwise search and update elements of best_results
		throughout.

		Parameters
		----------
		protocol: tuple of 2 ints
			First element is number of forward steps to take each iteration.
			The second is the number of reverse.  E.g., default is (2, 1).
			This results in two forward steps and one reverse step being taken
			each iteration.  E.g., if (1, 0) is passed, one forward step is
			taken followed by zero reverse steps.  That is, we do a forward
			search, etc.

		steps: int
			The total number of steps to take.  Note that if a step is
			requested but there are no mobile moves available no step will be
			taken.  This can happen, e.g., when requesting a forward step, but
			all mobile features are already in the predictor set, s.
		"""
		# write out full protocol for the current search request
		protocol = [i[0] for i in zip(
			cycle('f' * protocol[0] + 'r' * protocol[1]),
			range(steps))]

		# carry out stepwise search
		for p in protocol:
			if p == 'f':
				self._forward_step()
			else:
				self._reverse_step()
			self._update_best_results()

	def forward_cods(self):
		"""
		Returns the COD increase that would result from each possible movement
		of an element outside of the predictor set inside.

		Returns
		-------
		cod_gains : np.array (self.dimension, )
			This array's index i specifies the gain in target COD that would
			result if feature i were to move into s.  Values corresponding to
			unavailable moves are set to 0.
		"""
		cod_gains = np.full(self.dimension, 0, dtype=self.dtype)
		cod_gains[np.ix_(~self.s * self.mobile)] = np.einsum('ij,ij->i',
			(self.C - self.M)[np.ix_(
				~self.s * self.mobile, ~self.s * self.targets)],
			(self.C - self.M)[np.ix_(
				~self.s * self.mobile, ~self.s * self.targets)]) / (
				self.dtype(1.0) - self.C.diagonal()[np.ix_(
					~self.s * self.mobile)])
		return cod_gains

	def reverse_cods(self):
		"""
		Returns the COD decrease that would result from each possible movement
		of an element inside of the predictor set outside.

		Returns
		-------
		cod_costs : np.array (self.dimension, )
			This array's index i specifies the drop in target COD that would
			result if feature i were to move outside of s.  Values
			corresponding to unavailable moves are set to 0.
		"""
		cod_costs = np.full(self.dimension, 0, dtype=self.dtype)
		cod_costs[np.ix_(self.s * self.mobile)] = (
			self.targets[self.s * self.mobile]
			+ np.einsum('ij,ij->i',
				self.C[np.ix_(
					self.s * self.mobile, ~self.s * self.targets)],
				self.C[np.ix_(
					self.s * self.mobile, ~self.s * self.targets)])
			) / self.C.diagonal()[self.s * self.mobile]
		return cod_costs


def linFwd(bdtnp_file, dge_normalized_file):
	# bdtnp_file = "/Users/buutruong/TMTB/Computational_Research/UniSA/Dream_Challenges/DC_package/data/bdtnp.csv"
	# dge_normalized_file = "/Users/buutruong/TMTB/Computational_Research/UniSA/Dream_Challenges/DC_package/data/dge_normalized.csv"
	# insitu_df = pandas.read_csv(bdtnp_file)
	insitu_df = bdtnp_file
	marker_genes = list(insitu_df.columns)

	# Get expression of marker genes for cells
	# cell_df = pandas.read_csv(dge_normalized_file)
	cell_df = dge_normalized_file
	# all_genes = list(cell_df['Unnamed: 0'])  
	# cell_df = cell_df[cell_df.columns[1:]].T
	# cell_df.columns = all_genes
	
	cell_df = cell_df[marker_genes]
	# return cell_df.columns

	selector = FwdSelect()
	selector.fit(cell_df.values)
	idx = selector.ordered_features
	ls20 = [marker_genes[e] for e in idx][:20]
	ls40 = [marker_genes[e] for e in idx][:40]
	ls60 = [marker_genes[e] for e in idx][:60]
	return {'ls20':ls20, 'ls40':ls40, 'ls60':ls60}


def linRev(bdtnp_file, dge_normalized_file):
	insitu_df = bdtnp_file
	marker_genes = list(insitu_df.columns)

	# Get expression of marker genes for cells
	# cell_df = pandas.read_csv(dge_normalized_file)
	cell_df = dge_normalized_file
	# all_genes = list(cell_df['Unnamed: 0'])  
	# cell_df = cell_df[cell_df.columns[1:]].T
	# cell_df.columns = all_genes
	
	cell_df = cell_df[marker_genes]
	# return cell_df.columns

	selector = RevSelect()
	selector.fit(cell_df.values)
	idx = selector.ordered_features
	ls20 = [marker_genes[e] for e in idx][:20]
	ls40 = [marker_genes[e] for e in idx][:40]
	ls60 = [marker_genes[e] for e in idx][:60]
	return {'ls20':ls20, 'ls40':ls40, 'ls60':ls60}


def linGen(bdtnp_file, dge_normalized_file):
	# bdtnp_file = "./Data/bdtnp.csv"
	# dge_normalized_file = "./Data/dge_normalized1.csv"
	# insitu_df = pandas.read_csv(bdtnp_file)
	# return{"hello"}
	insitu_df = bdtnp_file
	marker_genes = list(insitu_df.columns)


	# Get expression of marker genes for cells
	# cell_df = pandas.read_csv(dge_normalized_file)
	cell_df = dge_normalized_file
	# all_genes = list(cell_df['Unnamed: 0'])  
	# cell_df = cell_df[cell_df.columns[1:]].T
	# cell_df.columns = all_genes
	
	cell_df = cell_df[marker_genes]
	# return cell_df.columns

	X = cell_df.values
	N = len(marker_genes)

	selector = GenSelect()
	selector.position(X)
	selector.search(protocol=(1, 0), steps=N)
	selector.search(protocol=(0, 1), steps=N)
	selector.search(protocol=(1, 0), steps=N)
	
	tmp = list(selector.best_results[20]['s']) 
	idx = [i for i in range(len(tmp)) if tmp[i]]
	ls20 = [marker_genes[e] for e in idx][:20]

	tmp = list(selector.best_results[40]['s']) 
	idx = [i for i in range(len(tmp)) if tmp[i]]
	ls40 = [marker_genes[e] for e in idx][:40]

	tmp = list(selector.best_results[60]['s']) 
	idx = [i for i in range(len(tmp)) if tmp[i]]
	ls60 = [marker_genes[e] for e in idx][:60]

	return {'ls20':ls20, 'ls40':ls40, 'ls60':ls60}

# def apb(a,b):
# 	return a+b
