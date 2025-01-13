function F=MyMatlabFunction1(f,Vars,var)
Vars=char(Vars);
if Vars(1)=='['
    Vars=Vars(2:end-1);
end
syms HELLO
Idx=has(f,var);
if length(Idx)>1
    f(Idx==0)=f(Idx==0)+HELLO;
elseif Idx==0
        f=f+HELLO;
end

f=char(f);
f=vectorize(f);
f=replace(f,'HELLO',strcat(['zeros(size(' char(var) '))']));
f=append(strcat(['@(' Vars ')']),f);
F=str2func(f);
end