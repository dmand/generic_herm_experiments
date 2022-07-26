
% Example 1: 
% Number of real eigenvalues of a family of random regular hermitian matrix
% pencils

clear all;

rng(5)
cc  = herm_gen_real(20,350);


plot(cc,'o','linewidth',1.5)
set(gca,'linewidth',1.5,'FontSize',16);
xlabel('Regular Hermitian matrix pencil') 
ylabel('Number of real eigenvalues')

%%

% Example 2: 
% Eigenvalues of singular hermitian matrix pencils

clear all;

n=17
r=9
number_of_ex = 1000

test = 1;

a = zeros(n,n);
b = zeros(n,n);

result_summary = zeros(2,floor(r/2)+1,number_of_ex);

for  ex = 1:number_of_ex
for s = 0:floor(r/2)
    for j = 1:s
        v = 1*rand(n,1)+1i*rand(n,1); 
        c = 1*rand(n,1)+1i*rand(n,1); 
        d = 1*rand(n,1)+1i*rand(n,1);
        a = a + v*c' + c*v';
        b = b + v*d' + d*v'; 
    end
    for k = 2*s+1:r
        u = 1*rand(n,1)+1i*rand(n,1); 
        a = a + 10*rand(1,1)*u*u';
        b = b + u*u';
    end
        eigenvalues = singgep(a,b); 
        if  norm(imag(eigenvalues)) > 1e-12*norm(a)*norm(b)
            %  testing if the eigenvalues are real
            test = 0
        end
        [n_of_eig,x]=size(eigenvalues);
        if ex == 1
        fprintf('s = %4.2f ; number of (real) eigenvalues =  %4.2f \n',s,n_of_eig)
        end
        result_summary(:,s+1,ex) = [s,n_of_eig]';
        a = zeros(n,n);
        b = zeros(n,n);
end
end

for jj = 1 : number_of_ex
    %  checking if in 10 000 experiments with the same s and l the numbers
    %  of eigenvalueswere the same
    if result_summary(:,:,ex) ~= result_summary(:,:,1)
        test  = 0;
    end
end
if test ==1
    fprintf('Same results were obtained for %4.2f randomly generated pencils \n',number_of_ex)
end 

%%

function [count] = herm_gen_real(m,n)
% Computes de eigenvalues of n random Hermitian mxm pencils and
% sums the number of real eigenvalues
count=zeros(n,1);
for s=1:n
    a = rherm(m)+log(s)^(1)*(s)*eye(m)/(n);
    b = rherm(m)+log(s)^(1)*(s)*eye(m)/(n);
    if isequal(a,a') && isequal(b,b') 
        e=eig(a,b);
            for j=1:m
             if abs(imag(e(j)))<=1e-10*(m*m*norm(norm(a),norm(b)))
                 count(s)=count(s)+1;
             end
            end
    else
       fprintf('Pencil is not Hermitian') 
    end
end
end

%%
function [ A ] = rherm(n)
% Generate a hermitian matrix with random values.
% rherm(n) returns a random hermitian matrix of size n*n
if exist('n')~=1, n=2; end
A = [zeros(1,0),(2*rand(1,1)-1),(2*rand(1,n-1)-1)+(2*rand(1,n-1)-1)*1i;zeros(n-1,n)];
for i=1:n-1, A = A + [zeros(n-i,n);zeros(1,n-i),(2*rand(1,1)-1),(2*rand(1,i-1)-1)+(2*rand(1,i-1)-1)*1i;zeros(i-1,n)];end
for i=1:n-1, for j=1:n-i, A = A + [zeros(n-i,n);zeros(1,j-1),conj(A(j,n-i+1)),zeros(1,n-j);zeros(i-1,n)];end;end
end
