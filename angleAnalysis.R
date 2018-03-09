# http://webspace.ship.edu/pgmarr/geo441/lectures/lec%2016%20-%20directional%20statistics.pdf


Angles = c(1, 359, 0, 341)*pi/180
mean(Angles)

Y = sin(Angles)
X = cos(Angles)

r = sqrt(mean(X)^2 + mean(Y)^2)

cab = mean(X)/r
sab = mean(Y)/r

meanAngle = atan(sab / cab)
meanAngle

