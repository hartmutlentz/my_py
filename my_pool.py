from multiprocessing import Pool

class ListPool(object):

	"""
	ListPool is a Pool of n Workers, where n is the number of processors

	ListPool needs

	liste = List of elements for which funktion should be calculated
	funktion = a function with two inputs (element, parameters)
			and two outputs (element, result_for_this_element)

	ListPool does

	run() = creates the Pool and evaluates function for every element in Liste
	result = returns a dictionary {element: result}
			result of course can be arbitrarly complicated

	ListPool cannot

	handle anything where the calculation of one element depends on the 
		result of another
	handle writing to common sources
	handle different function for different elements

	ListPool is not fully evaluated... it might have cause RAM-leakage...
	find out...
	"""

	def __init__(self, funktion, liste):
		#input
		self.liste = liste
		self.funktion = funktion
		#output
		self.__results = {}

	def cb(self, result):
		self.__results[result[0]] = result[1]

	def run(self):
		pool = Pool()
		for element in iter(self.liste):
			pool.apply_async(self.funktion,(element,), callback=self.cb)
		pool.close()
		pool.join()

	def getResults(self):
		return self.__results
	
	results = property(getResults)


