gene_number=3;
upper_limit=[0.5 5 50];
lower_limit=[0 0 0];
gene_size=[6 6 6]';
decimal=[2 2 2];
mutation_rate=0.01;
crossover_rate=0.7;
population_size=50;
generation=40;
loop=1;
g=1;
for i = 1:population_size
    for j = 1:gene_number
        gene(j,i,loop)=((rand*(upper_limit(j)-lower_limit(j)))+lower_limit(j));
    end
end
chromosome_size=sum(gene_size)+gene_number;
gene_start(1)=1;
for j=1:gene_number
gene_start(j+1)=gene_start(j)+gene_size(j)+1;
end
while loop <= generation
for i = 1:population_size
    for j=1:gene_number
        if gene(j,i,loop)>upper_limit(j)
            gene(j,i,loop)=upper_limit(j);
            elseif gene(j,i,loop) < lower_limit(j)
                gene(j,i,loop)=lower_limit(j);
        end
        if gene(j,i,loop) < 0
            pop(gene_start(j),i)=0;
        else
            pop(gene_start(j),i)=9;
        end
        temporary_gene(j,i)=abs(gene(j,i,loop));
        temporary_gene(j,i)=temporary_gene(j,i)/10^(decimal(j)-1);
        for k=gene_start(j):1:gene_start(j+1)-1
            pop(k,i)=temporary_gene(j,i)-rem(temporary_gene(j,i),1);
            temporary_gene(j,i)=(temporary_gene(j,i)-pop(k,i))*10;
        end
     end
end
a=1;
sumfitness=0;
for t=1:population_size
    m=0:0.01:10;
    numerator_1=0.752;
    denumerator_1=[0.157 1];
    sis1=tf(numerator_1,denumerator_1);
    numerator_2=[gene(1,a,loop) gene(2,a,loop) gene(3,a,loop)];
    a=a+1;
    denumerator_2=[1 0];
    sis2=tf(numerator_2,denumerator_2);
    sistop=series(sis1,sis2);
    new_sis=feedback(sistop,1);
    y=step(new_sis,m);
    error=0;
    for i=1:1001
        error=error+abs((y(i)-1)*0.01);
    end
objective_fitness(t)=error;
end
best_obj_function=objective_fitness(1);
for h=1:population_size
    if objective_fitness(h)<best_obj_function
        best_obj_function=objective_fitness(h);
    end
end
value(g)=best_obj_function;
g=g+1;
for t=1:population_size
    fitness(t)=1/objective_fitness(t);
    sumfitness = sumfitness + fitness(t);
end
[bestfitness(loop),bestmember]=max(fitness);
bestindividual(:,loop)=gene(:,bestmember,loop);
average_fitness(loop) = sumfitness/population_size;
for i=1:population_size
    pointer=rand*sumfitness
    members_number=1;
    toplam=fitness(1);
    while toplam < pointer
        members_number=members_number+1;
        toplam=toplam+fitness(members_number);
    end
    parent_chrom(:,i)=pop(:,members_number);
end
for s=1:population_size
    q=s;
    while q==s
        q = rand*population_size;
        q = q-rem(q,1)+1;
    end
    if crossover_rate > rand
        bolum = rand*chromosome_size;
        bolum = bolum-rem(bolum,1)+1;
        child(1:bolum,s)=parent_chrom(1:bolum,s);
        child(bolum+1:chromosome_size,s)=parent_chrom(bolum+1:chromosome_size,q);
    else
        child(:,s)=parent_chrom(:,s);
      end
end
for s=1:population_size
    for p=1:chromosome_size
        if mutation_rate > rand
            rand_gen=rand*10;
            while child(p,s)==rand_gen-rem(rand_gen,1)
                rand_gen=rand*10;
            end
            child(p,s)=rand_gen-rem(rand_gen,1);
        end
    end
end
pop=child;
loop=loop+1;
for i=1:population_size
    for j=1:gene_number
        gene(j,i,loop)=0;
        bookmark=1;
        for k=gene_start(j)+1:gene_start(j+1)-1
            place=decimal(j)-bookmark;
            gene(j,i,loop)=gene(j,i,loop)+(pop(k,i))*10^place;
            bookmark=bookmark+1;
        end
    end
end
end

t=0:0.01:10;
pay2=[response(1,1) response(2,1) response(3,1)];
payda2=[1 0];
sis2=tf(pay2,payda2);
sistop=series(sis1,sis2);
new_sis=feedback(sistop,1);
y=step(new_sis,t);
plot(t,y,'r');
hold on
pay_2=[0.029 2.84 37.6];
payda_2=[1 0];
sis_2=tf(pay_2,payda_2);
sis_top=series(sis1,sis2);
new_sis=feedback(sis_top,1);
[y,x]=step(new_sis_,t);
plot(t,y,'g-');