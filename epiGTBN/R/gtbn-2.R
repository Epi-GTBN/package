# Learn the structure of a Bayesian network using genetic tabu algorithm for 2-locus epistasis.

# Activate parameters including x, start, whitelist, blacklist, max.iter and debug
# Other parameters are deprecated and NOT RECOMMEND to use.

gtbn.2 = function(x, start, whitelist, blacklist, score, extra.args,
    restart, perturb, max.iter, maxp, optimized, debug = FALSE) {

    if(debug){
        cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ GTBN 2 NODES DEBUG MODE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
        cat("Activate parameters including x, start, whitelist, blacklist, max.iter and debug. Other parameters are deprecated and NOT RECOMMEND to use!\n")
		cat("\n============ initial start ============\n")
		print(start)
		cat("============ initial whitelist ============\n")
		print(whitelist)
		cat("============ initial blacklist ============\n")
		print(blacklist)
	}
	
	# cache nodes' labels.
	nodes = names(x)
	# cache the number of nodes.
	n.nodes = length(nodes)
	# need to avoid print() without enough output
	options(max.print=1000000)
	# set the iteration counter.
	iter = 0
	# check whether the score is score-equivalent.
	score.equivalence = is.score.equivalent(score, nodes, extra.args)
	# check whether the score is decomposable.
	score.decomposability = is.score.decomposable(score, extra.args)
	# allocate the cache matrix.
	cache = matrix(0, nrow = n.nodes, ncol = n.nodes)
	cache_ref = matrix(0, nrow = n.nodes, ncol = n.nodes)
	# nodes to be updated (all of them in the first iteration).
	updated = seq_len(n.nodes) - 1L

	# set the population size in all
    population.size = gtbn.population.size
	# the crossover probability 0.4~0.99
    cross.pro = gtbn.crossover.pro
	# mutation probability 0.0001~0.1
    mut.pro = gtbn.mut.pro
	# when the best individual not changed reached counter.unchanged times stop the ga 
	counter.unchanged = 4
	# the best individual has not changed for counter.unchanged.yet times
	counter.unchanged.yet = 0
	if(debug){
		cat("\n~~~~~~~~~~~~~~~~~~~~~ All the initiate value for GTBN are presented below ~~~~~~~~~~~~~~~~~~~~~\n")
            cat("\nnodes:", nodes,"\n");
            cat("n.nodes:", n.nodes,"\n");
            cat("score.equivalence:", score.equivalence,"\n");
            cat("score.decomposability:", score.decomposability,"\n");
            cat("updated:", updated,"\n");
			cat("population.size(default 50):", population.size,"\n");
			cat("cross.pro(default 0.7):", cross.pro,"\n");
			cat("mut.pro(default 0.002):", mut.pro,"\n");
			cat("counter.unchanged(default 4):", counter.unchanged,"\n");
			cat("counter.unchanged.yet(0):", counter.unchanged.yet,"\n");
        cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
	}
	
	# set the reference score.
	reference.score = per.node.score(network = start, score = score,
						targets = nodes, extra.args = extra.args, data = x)

	if(debug){
		cat("\n~~~~~~~~~~~~~~~~~~~~~ reference.score are presented below (network = start) ~~~~~~~~~~~~~~~~~~~~~\n")
			print(reference.score)
	}

	# convert the blacklist to an adjacency matrix for easy use.
	if (!is.null(blacklist))
		blmat = arcs2amat(blacklist, nodes)
	else
		blmat = matrix(0L, nrow = n.nodes, ncol = n.nodes)

	# convert the whitelist to an adjacency matrix for easy use.
	if (!is.null(whitelist))
		wlmat = arcs2amat(whitelist, nodes)
	else
		wlmat = matrix(0L, nrow = n.nodes, ncol = n.nodes)
	
	amat_ref = arcs2amat(start$arcs, nodes)
	.Call("score_cache_fill",
          nodes = nodes,
          data = x,
          network = start,
          score = score,
          extra = extra.args,
          reference = reference.score,
          equivalence = score.equivalence && optimized,
          decomposability = score.decomposability,
          updated = (if (optimized) updated else seq(length(nodes)) - 1L),
          amat = amat_ref,
		  cache = cache_ref,
          blmat = blmat,
          debug = debug)
	if(debug){
		cat("======== after score_cache_fill, cache_ref becomes ========\n")
			print(cache_ref)
	}
	# ---------------------beginning generate the primary individual-------------------------
    node1 = 0 
    node2 = 0 
    value = 0

    relation <- data.frame(node1, node2, value)

    myrbind.iter = 0 # instead of using rbind(), use data.frame() to do the same thing but faster.

    # transposition data matrix for binary use
    x_t = data.frame(t(x))
    x_classCol = x_t[nrow(x_t), ]
    x_t = x_t[-nrow(x_t), ]
    x_t = rbind(x_classCol, x_t)
    x_t = data.matrix(x_t)
    # rearrange colname order to put column `Class` in first place (stickler for now)
    nodes_new = c(nodes[length(nodes)], nodes[1 : (length(nodes)-1)])
    
    x_t.ncol = NCOL(x_t)
    x_t.nrow = NROW(x_t)
    .Call("num2binaryR",
            Rdata = x_t,
            RvarNum1 = 1, # column `Class` in 1st place
            Rncols = x_t.ncol,
            Rnrows = x_t.nrow)
    
    # time2start <- Sys.time()
    for(xi in 1 : (length(nodes_new)-1))
        for(yi in 1 : (length(nodes_new)-1))
            if(xi < yi) {
                myrbind.iter <- myrbind.iter + 1
                node1 <- c(nodes_new[xi+1])
                node2 <- c(nodes_new[yi+1])
                value <- mi.2(X = x_t, varNum2 = xi, varNum3 = yi)
                # cat("node1: ", node1, " node2: ", node2, " value: ", value, "\n")
                relation[myrbind.iter, 1] <- node1
                relation[myrbind.iter, 2] <- node2
                relation[myrbind.iter, 3] <- value
            }
    # time2end <- Sys.time()
    # time2 <- c(time2end - time2start)
    # cat("time full start: ", time2start, "time full end: ", time2end, "\n")
    # cat("time full iter: ", time2, "\n")
    
    # count after full iter before `repeat` time elapse
    # time4start <- Sys.time()
    
    # Sort by value
    o <- order(relation[, 3], decreasing = T)
    relation_order_all <- relation[o, ]
    # Remove the value column
    relation_order_all_noValCol <- (relation_order_all[, 1:2])
    # Take the largest 100
    relation_order_all_noValCol_core <- relation_order_all_noValCol[1:100,]
    
    if(debug){
        cat("\n======= top 100 relations are presented below ========\n")
        print(relation_order_all_noValCol_core)
    }

    # in case top 100 missed some nodes
    # 2 : n.nodes is to exclude 'Class' column which is row 1
	for (N_num in 2 : n.nodes) {
	  if ( sum(relation_order_all_noValCol_core == nodes_new[N_num]) == 0 ) {
	    # find the column labels of nodes that did not appear in the top 100
	    # then find the first line in the entire result and insert that interaction after the top 100 arcs
	    relation_order_all_noValCol_core = rbind(relation_order_all_noValCol_core, relation_order_all_noValCol[which(relation_order_all_noValCol == nodes_new[N_num])[1],])
	  }
	}

    if(debug){
        cat("\n======= top 100 plus missed nodes' relations, now presented below ========\n")
        print(relation_order_all_noValCol_core)
    }

    nodeInteration.amat = arcs2amat(as.matrix(relation_order_all_noValCol_core), nodes) 

    # chowliu.placeholder is to fix [Error in check.bn(x) : x must be an object of class 'bn'.].
    chowliu.placeholder = chow.liu(x) 


    # generate a acyclic network
    .Call("gnd_acyclic_network",
          amat1 = nodeInteration.amat,
          amat2 = amat_ref,
          nodes = nodes,
          wlmat = wlmat,
          blmat = blmat,
          debug = debug)

    if(debug){
        cat("\namat_ref is:\n")
            print(amat_ref)
    }
    amat(chowliu.placeholder) = amat_ref 
	
	if(debug){
		cat("\n~~~~~~~~~~~~~~~~~~~~~ score is presented below ~~~~~~~~~~~~~~~~~~~~~\n")
			print(score)
	}
    
    gtbn.start = chowliu.placeholder
    reference.score = per.node.score(network = gtbn.start, score = score,
						targets = nodes, extra.args = extra.args, data = x)
	if(debug){
		cat("\n~~~~~~~~~~~~~~~~~~~~~ reference.score are presented below ~~~~~~~~~~~~~~~~~~~~~\n")
			print(reference.score)

		cat("\n============================ network now is presented below ============================\n")
			print(gtbn.start)
	}
	# -----------------------------------------------------------------------
	if(debug){
        cat("\n============================ in score_cache_fill: ============================\n")
    }
    .Call("score_cache_fill",
          nodes = nodes,
          data = x,
          # network = start,
		  network = gtbn.start,
          score = score,
          extra = extra.args,
          reference = reference.score,
          equivalence = score.equivalence && optimized,
          decomposability = score.decomposability,
          updated = (if (optimized) updated else seq(length(nodes)) - 1L),
		  amat = amat_ref,
          cache = cache,
          blmat = blmat,
          debug = debug)
	if(debug){
		cat("\n======= after score_cache_fill, cache becomes ========\n")
			print(cache)
			# print(nodeInteration.amat)
	}

	if(debug)
	{
		cat("\n* the primary individual is:\n")
		print(gtbn.start)
	}
	
	# ---------------------ends of generate the primary individual-------------------------
	
	# the population used to storege all the individual 
	population = array(-1L, dim = c(n.nodes, n.nodes, population.size))
    
	
	# ---------------------beginning generate the primary polulation-------------------------
	.Call("set_rand_seed", debug = debug)
    for(i in 1 : population.size)
	{
		.Call("ga_generate_ind",
                   # amat = nodeInteration.amat,
                   amat = amat_ref, 
                   nodes = nodes,
                   wlmat = wlmat,
                   blmat = blmat,
                   debug = debug)
                   population[,,i] = amat_ref 
	}
	
	# ---------------------ends of generate the primary polulation-------------------------
	
	# init the next population
	population.next = population
	
	if(debug)
	{
		cat("\n* the primary population is:\n")
		print(population)
	}
	
	# keep the score of each individual
	population.score = array(-1, population.size)
	# used to calculate the score of each individual
	# use initial population to initialize network.temp
	network.temp = gtbn.start
	# record the best solution
	solution.best = -Inf
	# the cumpercentage of individual'score used in roulette selection
	individual.score.percentage = array(-1, population.size)
	# the cumpercentage of individual'score used in roulette selection
	individual.score.cumpercentage = array(-1, population.size)
	
	# set the length and allocate the tabu list for crossover.
	tabulist_length_crossover = gtbn.crossover.tabulist.length
	tabulist_crossover = vector("list", tabulist_length_crossover)
	current_crossover = 0
	# set the length and allocate the tabu list for mutation.
	tabulist_length_mutation = population.size
	tabulist_mutation = vector("list", tabulist_length_mutation)
	current_mutation = 0
	
    # time4end <- Sys.time()
    # time4 <- c(time4end - time4start)
    # cat("after full iter before repeat start: ", time4start, "after full iter before repeat end: ", time4end, "\n")
    # cat("after full iter before repeat: ", time4, "\n")


    # time3start <- Sys.time()
	repeat
	{
		if(debug)
		{
			cat("-------------",iter,"iteration-------------\n")
		}
		
		population = population.next
		
		
		# compute score of each individual in population
		for(i in 1 : population.size)
		{
			amat(network.temp) = population[,,i]
			population.score[i] = score(network.temp, x, type = score)
			if(debug)
				cat("compute score of each individual in population, right now is No.", i, " , score = "
					, population.score[i], "\n")
		}
		if(debug){
			cat("the score of each individual in population.\n")
			print(population.score)
		}
		
		# recored the score  of the best individual currently
		score.best = population.score[1]
		# recored the index of the best individual currently
		index.best = 1
		
		for(i in 2 : population.size)
		{
			if(population.score[i] > score.best)
				{
					index.best = i
					score.best = population.score[i]
					if(debug)
						cat("looks like when i = ", i, ", population.score[i]>score.best, update score.best\n")
				}
				if(debug)
					cat("No.", i, " in population.size(", population.size, "), score.best = ", score.best, "\n")
		}
		
		# put the best individual in the first place
		# exchange the highest score with the score of index 1
		score.temp = population.score[1]
		individual.temp = population[,,1]
		population.score[1] = population.score[index.best]
		population[,,1] = population[,,index.best]
		population.score[index.best] = score.temp
		population[,,index.best] = individual.temp
		
		if(debug){
			cat("the score of each individual in population after put the best in the first place.\n")
			print(population.score)
		}
		
        # one condition for stopping ga
		if(population.score[1] == solution.best)
		{
			# have not changed for a long time, confirm convergence
			if(counter.unchanged.yet >= counter.unchanged - 1)
			{
				if (debug)
					{
						cat("\n@ the best solution has not changed for a long time, stopping at iteration [ ")
						cat(counter.unchanged.yet, "]\n")
					}
					break;
					
			}
			counter.unchanged.yet = counter.unchanged.yet + 1
		}
		else
		{
		  # update score value
			solution.best = population.score[1]
			counter.unchanged.yet = 0
		}
		
		# another condition for stopping GTBN
		if(iter >= max.iter)
		{

			if (debug)
				{
					cat("\n@ reached max iteration, stopping at iteration: \n")
					cat(iter)
					cat("\n")
				}
			break;

		}
		
		# ---------------------beginning roulette selection-------------------------

		score.worst = min(population.score)*1.1 # make the worst possible to be choosen
		# each individual sub the worst score and then sum up
		score.all = sum(population.score) - population.size * score.worst
		if(debug){
			cat("score.worst",score.worst,"-------------------------\n")
			cat("score.all",score.all,"--------------------------\n")
		}
		if(score.all != 0)
		{
			# calculate percentage of each individual's score
			# for(i in 1 : population.size)
			# {
			# 	individual.score.percentage[i] = (population.score[i] - score.worst) / score.all
			# }
            individual.score.percentage = (population.score - score.worst) / score.all
			
            if(debug){
				cat("-------------------percentage------------------------------\n")
					print(individual.score.percentage)
			}
			percentage.temp = 0
			
            # calculate cumpercentage of each individual's score

			for(i in 1 : population.size)
			{
				percentage.temp = percentage.temp + individual.score.percentage[i]
				individual.score.cumpercentage[i] = percentage.temp
			}
			if(debug){		
				cat("-------------------individual.score.cumpercentage-------------\n")
					print(individual.score.cumpercentage)
			}
			population.next[,,1] = population[,,1]
			for(i in 2 : population.size)
			{
				pro = runif(1, 0, 1)
				for(j in 2 : population.size) # fixed the issue where this used to start with 1. It's 2 now.
				{
					if(individual.score.cumpercentage[j] >= pro && pro >= individual.score.cumpercentage[j-1])
					{
						if(debug)
							cat(i," times selection, ","pro:",pro," j:",j,"\n")
						population.next[,,i] = population[,,j]
						# cat(j)
						# cat("is selected\n")
						break
					}
				}
			}
			# ---------------------roulette selection ends-------------------------
			
		}
		else # The total score is 0, change to the next population
			polulation.next = population
		
		# ---------------------beginning crossover-------------------------
		# In this progress, this part is the only one that brings a cycle 
		
		for(i in 1 : population.size)
		{
			if(runif(1, 0, 1) < cross.pro)
			{
				# first individual choosed randomly to cross
				first.individual = round(runif(1, 1, population.size))
				# second individual choosed randomly to cross
				second.individual = round(runif(1, 1, population.size))
				
				if(first.individual == second.individual)
					next 
				
				# Calculate current item of tabu list for crossover
				current_crossover = as.integer(current_crossover %% tabulist_length_crossover)
				
				# ga_cross is in hc.cache.lookup.c
				.Call("ga_cross",
				      amat1 = population.next[,,first.individual],
				      amat2 = population.next[,,second.individual],
				      new_amat1 = population.next[,,first.individual],
				      new_amat2 = population.next[,,second.individual],
				      nodes = nodes,
				      wlmat = wlmat,
				      blmat = blmat,
				      tabu_list = tabulist_crossover,
				      current = current_crossover,
				      debug = debug)
				
				if(debug)
				{
					cat(first.individual, " individual has exchanged with ",
					second.individual, "\n")
				}
			}
		}
		# ---------------------crossover ends-------------------------

		
		# ---------------------beginning mutation-------------------------	
		# In this part, if a mutation bring a cycle, we refuse it, so the operation which bring a cycle is only in crossover 
		
		# Save all individuals in the current population into the tabu list
		for(i in 1 : population.size){
		  # update the value of current_mutation
		  current_mutation = as.integer((i-1) %% tabulist_length_mutation)
		  .Call("tabu_hash",
		        amat = population.next[,,i],
		        nodes = nodes,
		        tabu.list = tabulist_mutation,
		        current = current_mutation,
                debug = debug)
		}
		
		if(debug)
			cat("=================== beginning mutation ==================\n")
		for(i in 1 : population.size)
		{
			if(runif(1, 0, 1) < mut.pro)
			{
			  if(debug)
                  cat("stepping into mutation\n")
			  
			  # set up the score cache (BEWARE: in place modification!).
			  amat(network.temp) = population.next[,,i]
			  nparents = colSums(population.next[,,i])
			  if(debug) 
                  cat("caculate reference.score\n")
			  reference.score = per.node.score(network = network.temp, score = score,
			                                   targets = nodes, extra.args = extra.args, data = x)
			  
			  if(debug)
                  cat("call score_cache_fill\n")
			  .Call("score_cache_fill",
			        nodes = nodes,
			        data = x,
			        network = network.temp,
			        score = score,
			        extra = extra.args,
			        reference = reference.score,
			        equivalence = score.equivalence && optimized,
			        decomposability = score.decomposability,
			        updated = (if (optimized) updated else seq(length(nodes)) - 1L),
			        amat = population.next[,,i],
			        cache = cache,
			        blmat = blmat,
			        debug = debug)
			  
			  if(debug)
                  cat("to be added\n")
			  # select which arcs should be tested for inclusion in the graph (hybrid
			  # learning algorithms should hook the restrict phase here).
			  to.be.added = arcs.to.be.added(amat = population.next[,,i], nodes = nodes,
			                                 blacklist = blmat, whitelist = NULL, nparents = nparents,
			                                 maxp = maxp, arcs = FALSE)
			  
			  if(debug)
                  cat("call tabu_step")
			  # get the best arc addition/removal/reversal.
			  bestop = .Call("tabu_step",
			                 amat = population.next[,,i],
			                 nodes = nodes,
			                 added = to.be.added,
			                 cache = cache,
			                 reference = reference.score,
			                 wlmat = wlmat,
			                 blmat = blmat,
			                 tabu.list = tabulist_mutation,
			                 current = current_mutation,
			                 baseline = 0,
			                 nparents = nparents,
			                 maxp = maxp,
			                 debug = debug)
			  
			  
			  # find the index
			  mysearch.index <- function(nodes,nodes.length,target)
			  {
			    result = 0
			    for(i in 1 : nodes.length)
			    {
			      if(nodes[i] == target)
			      {
			        result = i
			        break
			      }
			    }
			    return (result)
			    
			  }
			  
			  if (bestop$op != FALSE) {
			    arcs.from = mysearch.index(nodes,n.nodes,bestop$from)
			    arcs.to = mysearch.index(nodes,n.nodes,bestop$to)
			    if(bestop$op == "set")
			    {
			      population.next[arcs.from,arcs.to,i] = 1L
			      if(debug > 0)
			      {
			        cat("for individual", i)
			        cat("adding", bestop$from, "->", bestop$to, ".\n")
			      }
			    }
			    else if(bestop$op == "drop")
			    {
			      population.next[arcs.from, arcs.to, i] = 0L 
			      if(debug > 0)
			      {
			        cat("for individual", i)
			        cat("droping", bestop$from, "->", bestop$to, ".\n")
			      }
			    }
			    else if(bestop$op == "reverse")
			    {
			      population.next[arcs.from,arcs.to,i] = 0L
			      population.next[arcs.to,arcs.from,i] = 1L
			      if(debug > 0)
			      {
			        cat("for individual", i)
			        cat("reversing", bestop$from, "->", bestop$to, ".\n")
			      }
			    }
			  }
			  
			  }
			}
		# ---------------------mutation ends-------------------------
		iter = iter + 1
	}
	result = gtbn.start
	amat(result) = population[,,1]
    
    # time3end <- Sys.time()
    # time3 <- c(time3end - time3start)
    # cat("time big repeat start: ", time3start, "time big repeat end: ", time3end, "\n")
    # cat("time big repeat: ", time3, "\n")

    if(!is.null(blacklist)){
        if(debug) cat("\ndealing with blacklist in final result...\n")
        for (blrow in 1 : nrow(blacklist)) {
            result = drop.arc(result, as.character(blacklist[blrow, 1]), as.character(blacklist[blrow, 2]))
        }
    }
    if(!is.null(whitelist)){
        if(debug) cat("\ndealing with whitelist in final result...\n")
        for (wlrow in 1 : nrow(whitelist)) {
            result = set.arc(result, as.character(whitelist[wlrow, 1]), as.character(whitelist[wlrow, 2]))
        }
    }
    
    if(debug){
        cat("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ final result is: ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
	        print(result)
        cat("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ $gy: bye! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
    }

	return (result)

}# GTBN.2NODES

