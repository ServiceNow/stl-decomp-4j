c data are 348 monthly values (29 years) in file co2
      real y(348), season(348), trend(348), rw(348), work(372,7)
      logical robust
      read(5,*)(y(i), i=1, 348)
      n=348
      np=12
      ns=35
      nt=19
      nl=13
      no=2
      ni=1
      nsjump=4
      ntjump=2
      nljump=2
      isdeg=1
      itdeg=1
      ildeg=1
      call stl(y,n,np,ns,nt,nl,isdeg,itdeg,ildeg,nsjump,ntjump,
     &nljump,ni,no,rw,season, trend, work)
       write(6,*) (season(i), i=1,348)
       write(6,*) (trend(i), i=1,348)
       write(6,*) (rw(i), i=1,348)
      robust=.true.
      call stlez(y,n,np,ns,isdeg,itdeg,robust,nconv,rw,season,trend,
     &work)
      write(6,*) nconv
       write(6,*) (season(i), i=1,348)
       write(6,*) (trend(i), i=1,348)
       write(6,*) (rw(i), i=1,348)
      end
