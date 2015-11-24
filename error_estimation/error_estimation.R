############################
##  likelihood: p(y|theta_1,theta_2)~N(theta1,1/theta2)
##  prior:      p(theta_1,theta_2) ~ p(theta_1|theta_2)p(theta_2)
##              p(theta_2)~gamma(a,b)
##              p(theta_1|theta_2)~N(mu,1/(tau*theta_2))
##              Then,
##              p(theta_1,theta_2)~normal-gamma(a,b,tau,mu)
##  posterior:  p(theta_1,theta_2|y)~normal-gamma
##
############################
