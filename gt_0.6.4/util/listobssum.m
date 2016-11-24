function listobssum(ts,te,rcvs,dirs,outf)
%
% list summary of observation data file
%
if nargin<1, ts=[2004,1,1]; end
if nargin<2, te=[2004,1,31]; end
if nargin<3, rcvs=prm_gpsrcvs; rcvs=rcvs(strcmp(rcvs(:,2),'IGS'),1); end
if nargin<4, dirs='k:\gps\obs_igs\%Y%m'; end
if nargin<5, outf='summary.txt'; end

td=caltomjd(ts):caltomjd(te);
s1{1}=sprintf('OBSERVATION DATA FILE SUMMARY');
s1{2}=sprintf('DATES     : %04d/%02d/%02d - %04d/%02d/%02d',ts,te);
s1{3}=sprintf('DIRECTORY : %s',dirs);
s1{4}=''; s1{5}='         '; s1{6}=s1{5}; s1{7}='----------';
for m=1:length(td)
    ep=mjdtocal(td(m));
    if ep(3)==1, s1{5}=[s1{5},' ']; s1{6}=[s1{6},' ']; s1{7}=[s1{7},'-']; end
    s1{5}=[s1{5},num2str(floor(ep(3)/10))];
    s1{6}=[s1{6},num2str(mod(ep(3),10))];
    s1{7}=[s1{7},'-'];
end
for n=1:length(rcvs)
    s2{n}=sprintf('%-6s : ',rcvs{n});
    rcv=rcvs{n}(max(1,end-3):end);
    for m=1:length(td)
        ep=mjdtocal(td(m)); doy=td(m)-caltomjd([ep(1),1,1])+1;
        if ep(3)==1, s2{n}=[s2{n},' ']; end
        dirr=strrep(dirs,'%Y',sprintf('%04d',ep(1)));
        dirr=strrep(dirr,'%m',sprintf('%02d',ep(2)));
        dirr=strrep(dirr,'%d',sprintf('%02d',ep(3)));
        dirr=strrep(dirr,'%s',rcv);
        file=fullfile(dirr,sprintf('%s%03d%0.%02d',rcv,doy,mod(ep(1),100)));
        fs=dir([file,'*']);
        if isempty(fs), s2{n}=[s2{n},'-']; else s2{n}=[s2{n},'o']; end
    end
end
f=fopen(outf,'wt'); if f<0, return, end
for n=1:length(s1), fprintf(f,'%s\n',s1{n}); end
for n=1:length(s2), fprintf(f,'%s\n',s2{n}); end
fclose(f);
