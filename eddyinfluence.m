dp=0.5/10;delta=3.2;m=10000;n=500;simu=25;
aa=1;
for i=1:simu
    Y(:,aa:m*i)=YY{i}(1:n,:);
    aa=m*i+1;
end
aa=1;
for i=1:simu
    meet(:,aa:m*i)=meett{i}(1:n,:);
    aa=m*i+1;
end
aa=1;
for i=1:simu
    drx(:,aa:m*i)=ddrx{i}(1:n,:);
    aa=m*i+1;
end
aa=1;
for i=1:simu
    dry(:,aa:m*i)=ddry{i}(1:n,:);
    aa=m*i+1;
end

m=1;n=1;
for j=1:10000
    for i=1:500
        if meet(i,j)==1
            drxa(m)=drx(i,j);
            drya(m)=dry(i,j);
            ya(m)=Y(i,j);
            m=m+1;
        elseif meet(i,j)==2
            drxb(n)=drx(i,j);
            dryb(n)=dry(i,j);
            yb(n)=Y(i,j);
            n=n+1;
        end
    end
end

k=50;
c_y=linspace(dp,delta,k);
for i=1:k-1
    clear a
    a(1,:)=find(ya>c_y(i)&ya<c_y(i+1));
    m=1;
    for j=1:length(a)
        drxaa(m)=drxa(a(j));
        dryaa(m)=drya(a(j));
        m=m+1;
    end
    clear b
    b(1,:)=find(yb>c_y(i)&yb<c_y(i+1));
    n=1;
    for j=1:length(b)
        drxbb(n)=drxb(b(j));
        drybb(n)=dryb(b(j));
        n=n+1;
    end
    vxa(i,1)=var(drxaa,0,2);
    vya(i,1)=var(dryaa,0,2);
    vxb(i,1)=var(drxbb,0,2);
    vyb(i,1)=var(drybb,0,2);
    y(i,1)=(c_y(i)+c_y(i+1))/2;
    clear drxaa dryaa drxbb drybb
end
figure(1);plot(vxa,y,'b.','markersize',20);hold on;plot(vxb,y,'r.','markersize',20);
title('streamwise variance of the turbulence term');xlabel('variance');ylabel('Y(cm)');legend('type-A eddies','type-B eddies')
figure(2);plot(vya,y,'b.','markersize',20);hold on;plot(vyb,y,'r.','markersize',20);
title('vertical variance of the turbulence term');xlabel('variance');ylabel('Y(cm)');legend('type-A eddies','type-B eddies')


