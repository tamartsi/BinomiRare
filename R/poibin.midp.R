poibin.midp <-
function(n.carrier, n.D.carrier, prob.vec){
	stopifnot(n.D.carrier <= n.carrier, length(prob.vec) == n.carrier)
	d.poibin <- dpoibin(0:n.carrier, prob.vec)
	prob.cur <- d.poibin[n.D.carrier + 1]
	mid.p <- 0.5*prob.cur + sum(d.poibin[d.poibin < prob.cur])
	return(mid.p)
}
