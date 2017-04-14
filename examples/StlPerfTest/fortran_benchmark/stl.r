subroutine stl(y,n,np,ns,nt,nl,isdeg,itdeg,ildeg,nsjump,ntjump,nljump,ni,no,rw,
	season,trend,work)

integer n, np, ns, nt, nl, isdeg, itdeg, ildeg, nsjump, ntjump, nljump, ni, no, k
integer newns, newnt, newnl, newnp
real y(n), rw(n), season(n), trend(n), work(n+2*np,5)
logical userw

userw = .false.
k = 0
do i = 1,n
	trend(i) = 0.0
newns = max0(3,ns)
newnt = max0(3,nt)
newnl = max0(3,nl)
newnp = max0(2,np)
if(mod(newns,2)==0) newns = newns + 1	#make odd
if(mod(newnt,2)==0) newnt = newnt + 1
if(mod(newnl,2)==0) newnl = newnl + 1
repeat {
	call onestp(y,n,newnp,newns,newnt,newnl,isdeg,itdeg,ildeg,nsjump,ntjump,nljump,
	ni,userw,rw,season, trend, work)
	k = k+1
	if(k>no) break
	do i = 1,n
		work(i,1) = trend(i)+season(i)
	call rwts(y,n,work(1,1),rw)
	userw = .true.
	}
if (no<=0)
	do i = 1,n
		rw(i) = 1.0
return
end



subroutine ess(y,n,len,ideg,njump,userw,rw,ys,res)

integer n, len, ideg, njump, newnj, nleft, nright, nsh, k, i, j
real y(n), rw(n), ys(n), res(n), delta
logical ok, userw

if(n<2) { ys(1) = y(1); return }
newnj = min0(njump, n-1)
if(len>=n) {	#len > or = n
	nleft = 1
	nright = n
	do i = 1,n,newnj {
		call est(y,n,len,ideg,float(i),ys(i),nleft,nright,res,userw,rw,ok)
		if(!ok) ys(i) = y(i)
		}
	}
else if(newnj==1) {	# newnj equal to one, len less than n
	nsh = (len+1)/2
	nleft = 1
	nright = len
	do i = 1,n {	# fitted value at i
		if(i>nsh && nright!=n) {
			nleft = nleft+1
			nright = nright+1
			}
		call est(y,n,len,ideg,float(i),ys(i),nleft,nright,res,userw,rw,ok)
		if(!ok) ys(i) = y(i)
		}
	}
else { # newnj greater than one, len less than n
	nsh = (len+1)/2
	do i = 1,n,newnj {	# fitted value at i
		if(i<nsh) {
			nleft = 1
			nright = len
			}
		else if(i>=n-nsh+1) {
			nleft = n-len+1
			nright = n
			}
		else {
			nleft = i-nsh+1
			nright = len+i-nsh
			}
		call est(y,n,len,ideg,float(i),ys(i),nleft,nright,res,userw,rw,ok)
		if(!ok) ys(i) = y(i)
		}
	}
if(newnj!=1){
	do i = 1,n-newnj,newnj {
		delta = (ys(i+newnj)-ys(i))/float(newnj)
		do j = i+1,i+newnj-1
			ys(j) = ys(i)+delta*float(j-i)
		}
	k = ((n-1)/newnj)*newnj+1
	if(k!=n) {
		call est(y,n,len,ideg,float(n),ys(n),nleft,nright,res,userw,rw,ok)
		if(!ok) ys(n) = y(n)
		if(k!=n-1) {
			delta = (ys(n)-ys(k))/float(n-k)
			do j = k+1,n-1
				ys(j) = ys(k)+delta*float(j-k)
			}
		}
	}
return
end


subroutine est(y,n,len,ideg,xs,ys,nleft,nright,w,userw,rw,ok)

integer n, len, ideg, nleft, nright, j
real y(n), w(n), rw(n), xs, ys, range, h, h1, h9, a, b, c, r
logical userw,ok

range = float(n)-float(1)
h = amax1(xs-float(nleft),float(nright)-xs)
if (len>n) h = h+float((len-n)/2)
h9 = .999*h
h1 = .001*h
# compute weights
a = 0.0
do j = nleft,nright {
	w(j) = 0.
	r = abs(float(j)-xs)
	if (r<=h9) {
		if (r<=h1) w(j) = 1.
		else w(j) = (1.0-(r/h)**3)**3
		if (userw) w(j) = rw(j)*w(j)
		a = a+w(j)
		}
	}
if (a<=0.0)
	ok = .false.
else {	# weighted least squares
	ok = .true.
	do j = nleft,nright	# make sum of w(j) == 1
		w(j) = w(j)/a
	if ((h>0.)&(ideg>0)) {	# use linear fit
		a = 0.0
		do j = nleft,nright	# weighted center of x values
			a = a+w(j)*float(j)
		b = xs-a
		c = 0.0
		do j = nleft,nright
			c = c+w(j)*(float(j)-a)**2
		if (sqrt(c)>.001*range) {
			b = b/c
# points are spread out enough to compute slope
			do j = nleft,nright
				w(j) = w(j)*(b*(float(j)-a)+1.0)
			}
		}
	ys = 0.0
	do j = nleft,nright
		ys = ys+w(j)*y(j)
	}
return
end



subroutine fts(x,n,np,trend,work)

integer n, np
real x(n), trend(n), work(n)

call ma(x,n,np,trend)
call ma(trend,n-np+1,np,work)
call ma(work,n-2*np+2,3,trend)
return
end


subroutine ma(x, n, len, ave)

integer n, len, i, j, k, m, newn
real x(n), ave(n), flen, v

newn = n-len+1
flen = float(len)
v = 0.0
# get the first average
do i = 1,len
	v = v+x(i)
ave(1) = v/flen	
if (newn>1) {
	k = len
	m = 0
	do j = 2, newn {
# window down the array
		k = k+1
		m = m+1
		v = v-x(m)+x(k)
		ave(j) = v/flen	
	}
}
return
end


subroutine onestp(y,n,np,ns,nt,nl,isdeg,itdeg,ildeg,nsjump,ntjump,nljump,ni,
	userw,rw,season,trend,work)

integer n,ni,np,ns,nt,nsjump,ntjump,nl,nljump,isdeg,itdeg,ildeg
real y(n),rw(n),season(n),trend(n),work(n+2*np,5)
logical userw

do j = 1,ni {
	do i = 1,n
		work(i,1) = y(i)-trend(i)
	call ss(work(1,1),n,np,ns,isdeg,nsjump,userw,rw,work(1,2),work(1,3),work(1,4),
		work(1,5),season)
	call fts(work(1,2),n+2*np,np,work(1,3),work(1,1))
	call ess(work(1,3),n,nl,ildeg,nljump,.false.,work(1,4),work(1,1),work(1,5))
	do i = 1,n
		season(i) = work(np+i,2)-work(i,1)
	do i = 1,n
		work(i,1) = y(i)-season(i)
	call ess(work(1,1),n,nt,itdeg,ntjump,userw,rw,trend,work(1,3))
	}
return
end


subroutine rwts(y,n,fit,rw)

integer mid(2), n
real y(n), fit(n), rw(n), cmad, c9, c1, r

do i = 1,n
	rw(i) = abs(y(i)-fit(i))
mid(1) = n/2+1
mid(2) = n-mid(1)+1
call psort(rw,n,mid,2)
cmad = 3.0*(rw(mid(1))+rw(mid(2)))	#6 * median abs resid
c9 = .999*cmad
c1 = .001*cmad
do i = 1,n {
	r = abs(y(i)-fit(i))
	if (r<=c1) rw(i) = 1.
	else if (r<=c9) rw(i) = (1.0-(r/cmad)**2)**2
	else rw(i) = 0.
	}
return
end



subroutine ss(y,n,np,ns,isdeg,nsjump,userw,rw,season,work1,work2,work3,work4)

integer n, np, ns, isdeg, nsjump, nright, nleft, i, j, k
real y(n), rw(n), season(n+2*np), work1(n), work2(n), work3(n), work4(n), xs
logical userw,ok

for(j=1; j<=np; j=j+1){
	k = (n-j)/np+1
	do i = 1,k
		work1(i) = y((i-1)*np+j)
	if(userw)
		do i = 1,k
			work3(i) = rw((i-1)*np+j)
	call ess(work1,k,ns,isdeg,nsjump,userw,work3,work2(2),work4)
	xs = 0
	nright = min0(ns,k)
	call est(work1,k,ns,isdeg,xs,work2(1),1,nright,work4,userw,work3,ok)
	if(!ok) work2(1) = work2(2)
	xs = k+1
	nleft = max0(1,k-ns+1)
	call est(work1,k,ns,isdeg,xs,work2(k+2),nleft,k,work4,userw,work3,ok)
	if(!ok) work2(k+2) = work2(k+1)
	do m = 1,k+2
		season((m-1)*np+j) = work2(m)
	}
return
end
subroutine stlez(y, n, np, ns, isdeg, itdeg, robust, no, rw, season, trend, work)

logical robust
integer n, i, j, np, ns, no, nt, nl, ni, nsjump, ntjump, nljump, newns, newnp
integer isdeg, itdeg, ildeg
real y(n), rw(n), season(n), trend(n), work(n+2*np,7)
real maxs, mins, maxt, mint, maxds, maxdt, difs, dift

ildeg = itdeg
newns = max0(3,ns)
if(mod(newns,2)==0) newns = newns+1
newnp = max0(2,np)
nt = (1.5*newnp)/(1 - 1.5/newns) + 0.5
nt = max0(3,nt)
if(mod(nt,2)==0) nt = nt+1
nl = newnp
if(mod(nl,2)==0) nl = nl+1
if(robust) ni = 1
else ni = 2
nsjump = max0(1,int(float(newns)/10 + 0.9))
ntjump = max0(1,int(float(nt)/10 + 0.9))
nljump = max0(1,int(float(nl)/10 + 0.9))
do i = 1,n
	trend(i) = 0.0
call onestp(y,n,newnp,newns,nt,nl,isdeg,itdeg,ildeg,nsjump,ntjump,nljump,ni,
	.false.,rw,season,trend,work)
no = 0
if(robust){
	for(j=1; j<=15; j=j+1){	#robustness iterations
		do i = 1,n{	#initialize for testing
			work(i,6) = season(i)
			work(i,7) = trend(i)
			work(i,1) = trend(i)+season(i)
			}
		call rwts(y,n,work(1,1),rw)
		call onestp(y, n, newnp, newns, nt, nl, isdeg, itdeg, ildeg, nsjump,
			ntjump, nljump, ni, .true., rw, season, trend, work)
		no = no+1
		maxs = work(1,6)
		mins = work(1,6)
		maxt = work(1,7)
		mint = work(1,7)
		maxds = abs(work(1,6) - season(1))
		maxdt = abs(work(1,7) - trend(1))
		do i = 2,n{
			if(maxs<work(i,6)) maxs = work(i,6)
			if(maxt<work(i,7)) maxt = work(i,7)
			if(mins>work(i,6)) mins = work(i,6)
			if(mint>work(i,7)) mint = work(i,7)
			difs = abs(work(i,6) - season(i))
			dift = abs(work(i,7) - trend(i))
			if(maxds<difs) maxds = difs
			if(maxdt<dift) maxdt = dift
			}
		if((maxds/(maxs-mins)<.01) & (maxdt/(maxt-mint)<.01)) break
		}
	}
if(!robust){
	do i = 1,n
		rw(i) = 1.0
	}
return
end
subroutine psort(a,n,ind,ni)
real a(n)
integer n,ind(ni),ni
integer indu(16),indl(16),iu(16),il(16),p,jl,ju,i,j,m,k,ij,l
real t,tt
if (n<0||ni<0) return #ERROR(number of items<0)
if(n<2 | ni==0) return
jl = 1
ju = ni
indl(1) = 1
indu(1) = ni
#  arrays indl, indu keep account of the portion of ind related to the
#  current segment of data being ordered.
i = 1
j = n
m = 1
repeat {
	if (i<j) go to 10
	repeat {
		m = m-1
		if (m==0) break 2
		i = il(m)
		j = iu(m)
		jl = indl(m)
		ju = indu(m)
		if (jl<=ju) {
			while (j-i>10) {
#  first order a(i),a(j),a((i+j)/2), and use median to split the data
				10  k = i
				ij = (i+j)/2
				t = a(ij)
				if (a(i)>t) {
					a(ij) = a(i)
					a(i) = t
					t = a(ij)
					}
				l = j
				if (a(j)<t) {
					a(ij) = a(j)
					a(j) = t
					t = a(ij)
					if (a(i)>t) {
						a(ij) = a(i)
						a(i) = t
						t = a(ij)
						}
					}
				repeat {
					l = l-1
					if (a(l)<=t) {
						tt = a(l)
						repeat
#  split the data into a(i to l).lt.t, a(k to j).gt.t
							k = k+1
							until(a(k)>=t)
						if (k>l) break
						a(l) = a(k)
						a(k) = tt
						}
					}
				indl(m) = jl
				indu(m) = ju
				p = m
				m = m+1
#  split the larger of the segments
				if (l-i<=j-k) {
					il(p) = k
					iu(p) = j
					j = l
					repeat {
						if (jl>ju) next 3
						if (ind(ju)<=j) break
						ju = ju-1
						}
					indl(p) = ju+1
					}
				else {
					il(p) = i
					iu(p) = l
					i = k
					repeat {
#  skip all segments not corresponding to an entry in ind
						if (jl>ju) next 3
						if (ind(jl)>=i) break
						jl = jl+1
						}
					indu(p) = jl-1
					}
				}
			if (i==1) break
			i = i-1
			repeat {
				i = i+1
				if (i==j) break
				t = a(i+1)
				if (a(i)>t) {
					k = i
					repeat {
						a(k+1) = a(k)
						k = k-1
						}
						until(t>=a(k))
					a(k+1) = t
					}
				}
			}
		}
	}
return
end
