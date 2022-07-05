% leastsqTutorial.m  script
%
% 
% This is a beginner's tutorial.  If you play with the concepts in this 
% tutorial, you can go from knowing next to nothing about linear algebra 
% to an intuition that will help with multiple linear regression.
%
% M.N. Shadlen for Cold Spring Harbor 1998
% return

clear all; close all;

% suppose you have ordered pairs
data = [-1 1 2; 1 1 3]'

% Our job is to find the best fitting line for these three points.
% Here are the points
hf1 = figure(1), clf
plot(data(:,1), data(:,2),'o')
set(gca,'XLim',[-3 3],'YLim',[-1 4])


% We are going to set this up as a matrix equation.
% Ax = b.  The term x is somewhat unfortunate because x is 
% not the independent variable, but the constant and slope. 
% b is almost as unfortunate since many of us think of b as the
% slope of a line.  Here it is the value plotted on the ordinate.
% Let's start by making A.  It is a 3 x 2  matrix 
% of the left side variables p. 161

A = [ones(3,1) data(:,1)]

% Notice that the first column of A is just a column of ones. 
% The second column contains the values plotted on the abscissa -- what
% we would normally call the x-variable: the values -1 1 2.
% To avoid confusion, let's call the "x-values" a

a = data(:,1)

% let's call the "y-values" vector b
b = data(:,2)

% We know A and we know b.  We want to find the best solution for x. 
% The first thing you need to convince yourself of is that the 2 by 1 column 
% vector, x, contains the constant and slope of the best fitting line.
% This is a good time to remind yourself of how a matrix multiplies a vector.

% Let's take a wild guess at x.  Let's guess that the intercept is .5 and
% the slope is 1
x = [.5 1]'

% Notice that you get out 3 values from A*x
b_guess = A*x
hold on, plot(a,b_guess,'+',a,b_guess,':')

% That's a pretty crummy fit, but that's not the point.  All you're supposed
% to see right now is that by going across the rows of A and down the
% column vector x, you get the equation 
[	A(1,1)*x(1) + A(1,2)*x(2)
	A(2,1)*x(1) + A(2,2)*x(2)
	A(3,1)*x(1) + A(3,2)*x(2)
]

% If you don't understand that last equation, stop.  There's no point
% in going any further. Ask for help or pick up a book.  
% If you understand the equation, then ask yourself the following.
% Suppose we made A as a column of 1s a column like A(:,2) and a third
% column of the squares:  A(:,3) = A(:,2).*A(:,2).  We would be looking for 
% a column vector, x, that has 3 components.  Do you see that these are the
% coefficients of the best fitting quadratic? x1 + x2*a + x3*a^2 = b
% Best to convince yourself of that too before going on.

% Now this way of setting up the equation jives with your sense of 
% matrix times vector by going along the rows of the matrix and multiplying
% the elements by those going down the column vector.  Another way to put
% that is that each element in b is a dot-product of the corresponding
% row in A with the column vector x. 
% That's worth thinking about, but it does not lend itself to a particularly
% lucid view of the regression problem. 
% For this, it is worth thinking about the matrix times the vector in
% a different way.  Instead of going along the rows of A, consider
% the columns.  This is actually the easier way to think about it.
% Notice that b is the same dimension as a column in A (3 by 1).
% Convince yourself that b is just the weighted sum of the columns in A.
% The weights come from the elements of x.

b_guess = x(1)*A(:,1) + x(2)*A(:,2)


% Again, it's not worth going on if you don't see this version of A*x

% Now here is where things start getting interesting.
% Once we recognize that b is the weighted sum of just 2 vectors, 
% it must be the case that b lies in the plane that these vectors 
% define.  No matter how you cut it, you can only find vectors, b, that you 
% can get to by adding the two vectors that form the columns of A.
% This plane is called the column-space of A.  It includes the origin.

% It helps to look at a picture.  The way I'm going to draw this picture
% includes a step or two that we're not ready for.  Just execute it
% for the time being.  Later, we'll come back to why the matrix algebra
% works.


% Let's start by drawing our two column vectors
hf2 = figure(2)
hold off; hp = plot3(A(1,:),A(2,:),A(3,:),'ro','MarkerFaceColor','r'), hold on;
hl1= line([0 A(1,1)], [0 A(2,1)], [0 A(3,1)], 'Color', 'r')
hl2=line([0 A(1,2)], [0 A(2,2)], [0 A(3,2)],'Color', 'r')
s = sprintf(' text(%d,%d,%d,''(%d,%d,%d)'')', [A(:,1) A(:,1)]), eval(s)
s = sprintf(' text(%d,%d,%d,''(%d,%d,%d)'')', A(:,[2 2])), eval(s)
set(gca,'Box', 'on','XLim',[-4 4],'YLim',[-4 4],'ZLim',[-2 3])
set(gca,'XGrid','on','YGrid','on','ZGrid','off')

% You can imagine the plane
hfill = fill3([0 A(1,:) 0], [0 A(2,:) 0], [0 A(3,:) 0], 'b','EraseMode','xor')

% Now here's the part that you need to ignore for a moment.
% The next line makes a projection matrix. We're going to use it 
% to project the xy plane into the column space of A.
P = A * inv(A'*A) * A'
q = P * 2 * [1 1 1; 1 -2 1; -2 -2 1; -2 1 1; 1 1 1]'
h2 = fill3(q(1,:), q(2,:), q(3,:) , 'b','EraseMode','xor','EdgeColor','r')
delete(hfill)
xlabel('b1'), ylabel('b2'), zlabel('b3')

figure(2)
for k = 332:-4:300
    set(gca,'View',[k 30])
			drawnow
end
for k = 30:5:50
  set(gca,'View',[300 k])
	drawnow
end

% Are you convinced that no matter what we guess for our solution, x, we
% will end up with something in the plane?  You should have convinced yourself
% of this earlier.  But here is a quick demo.


% It is important to realize that any guess we make for the answer to
% our problem will make a vector in the plane spanned by A's column vectors.
% We'll make 5 random guesses for g and see what we get out
for i = 1:5
	g = 2 * (rand(2,1) - .5)
	res = A * g
	hg(i) = plot3(res(1),res(2),res(3),'c.')
	hl(i) = line([0 res(1)], [0 res(2)], [0 res(3)], 'Color','y')
end

% You might want to replay the axis-rotation again. Depending on where
% those random vectors point, it may or may not be obvious that they are in 
% the plane.

% Now that you have convinced yourself that no matter what we
% choose for x, we are going to get out (after multiplication by A) a 
% vector that lies in the plane. 
% That is unfortunate because the perfect solution would be the
% vector b, and it is usually the case that b does not lie in our plane!
% What a bummer.  But then again, if it did, we would have a perfect
% solution -- and that's not what we're dealing with in least squares.

% Let's look at the vector b.  Remember this is the list of ordinate
% values for the points in Figure 1.

plot3(b(1),b(2),b(3),'g*')
line([0 b(1)], [0 b(2)], [0 b(3)], 'Color','g')
set(gca,'ZLim',[-1 3])

% It helps to get rid of the extra lines
delete(hl)
delete(hg)

% It also helps to anchor the orgin 

line([0 0],[0 0],[-1 0],'LineWidth',3,'Color','r')

% and a little spinning couldn't hurt

figure(2)
k = 50
i = 300
set(gca,'View',[i k])
drawnow, % pause
for k = 50:-4:30
 	set(gca,'View',[i k])
	[i k]
	drawnow
	% pause
end
for i = 300:-5:220
 	set(gca,'View',[i k])
	[i k]
	drawnow
	% pause
end

% So now what?  Clearly b is not in the column space of A.
% Yet no matter what we choose for x, we will get a set
% of predicted values that form a vector that is in the column-space.
% Of course we want to find x so that the vector we get out is the
% closest one to the real solution.  
% Let's call the solution, x_bar.  After all, there is no vector x, such
% that A*x = b.

% Now from the picture in front of you, it ought to be obvious that the
% closest we can get to b is a vector that ends in the projection of b on
% the blue plane.  I'm going to drop the end of vector b onto the plane.
% Try to ignore the math for the moment and grasp the geometry.  

bproj = A * (inv(A' * A) * A' * b)
hp = line([0 bproj(1)], [0 bproj(2)], [0 bproj(3)])
set(hp,'Color','g','LineStyle','--','LineWidth',3)
% delete(hp)

% The thick broken line is the projection of b onto the column-space of
% A.  Let's call this vector bproj.  It is the predicted value that we would
% get out of A*x_bar = bproj.  One way to see this is to return to the
% first plot and look at where bproj lies. These are the green asterisks 
% connected by the green broken line.  

figure(1)
hold off, plot(a,b,'o'), hold on
plot(a,bproj,'g*', a,bproj,'g--')
set(gca,'XLim',[-3 3],'YLim',[-1 4])

% Our bproj is the column vector containing the ordinate values of these 
% three points.  You can see that these are the predicted values that lie
% on the best fitting line.  The way that we would say that this is the best 
% fitting line is that the distance from the asterisks to the data values 
% (vertically) is minimized.

% This distance is the difference between b and bproj.
% Let's return to the vector diagram and convince ourselves.

% Visualizing this error vector is the key to understanding least squares.

figure(2)

% It helps to construct a line that represents the error between b and
% bproj.  The error, E, is

E = b - bproj

% Rather than displaying this vector emanating from the origin, we
% will add it to the end of bproj

hE = line([bproj(1) b(1)], [bproj(2) b(2)], [bproj(3) b(3)])
set(hE,'Color','y','LineStyle','--','LineWidth',3)

% Do you appreciate that this E vector is orthogonal to the plane?
% There are lots of ways to talk about this, but you have to see
% and believe it first.  Don't go on unless you are absolutely comfortable
% with the idea that E is orthogonal to the plane spanned by the columns of 
% A.

% O.K., you made it.  Do you see that regardless of the dimension of the 
% problem, the least squares solution will be a projection of some
% high dimensional vector, b, on some lower dimensional column space of A.  
% Sticking with the 2-dimensional case, if we were fitting a
% line to 4 points instead of 3, we would be projecting a 4-dimensional 
% vector, b, into the 2-dimensional plane spanned by the 2 4-dimensional
% column vectors of A.  One of these would be [1 1 1 1]'  

% Let's return to the matter of E and how to talk about it.  
% We have already said that it is orthogonal to the columns of A.
% That means that the dot product of E with any column of A should be 
% zero.  Since E is a column vector, we either have to take its transpose 
% and multiply E'*A(:,1) 

E' * A(:,1)     % 1st column 
E' * A(:,2)			% 2nd column

% or take the transpose of the column vectors and multiply these times
% the E column vector.

A(:,1)' * E		  % 1st column 
A(:,2)' * E			% 2nd column

% I'm assuming you know about dot products.  If you don't, you'll have
% to pick up a book. You might just recall that it is the product 
% of the lengths times the cosine of the angle between the vectors.
% So z'*z is the square of the length of z.  You can easily verify
% that the dot product of two orthogonal vectors is 0
% e.g., visualize and compute

[1 0] * [0 1]'
[1 1] * [-1 1]'
[1 1] * [-2.7 2.7]'

% ... and so forth.
 
% If you don't understand dot products (sometimes called inner products)
% don't go on. 

% If you understand the notion of inner products and the fact that E is
% orthogonal to the columnspace of A, then we can move on.
% The first thing to say is that we can combine the dot products above
% into a single matrix operation.  Starting with E transpose, 

E' * A		% should be a row vector with 2 0's. (You get very tiny numbers)

% Or we can transpose A and multiply it times the E column vector

A' * E		% makes a column vector with 2 zeros. (within floating point precision)


% There is a very imporant idea here.  The Error vector, which is, after all,
% the residuals between the observed data and the fit, is orthogonal to
% the columnspace of A.  It lies in the left nullspace of A.  The reason 
% this is called the left nullspace is because we have to left-multiply 
% E * A to get our zeros.

% Once we see that this is so, we can derive the NORMAL EQUATIONS
% in matrix form.  We are looking for the best approximation to
% x in the equation A*x = b.  We have already admitted that there 
% is no vector x that will work because A*x always lies in the columnspace
% of A, whereas b does not.  So, instead, we'll solve for xbar 
% We want the best solution to A*xbar = b.   Again, we admit we can't 
% achieve equality.  So we say that what we really want is to minimize the
% error between our fit, A*xbar, and the observed values, b.  
% In least squares, we want to minimize the  sum of the (b - A*xbar).^2
% What you should be absolutely convinced of is that 
% b - A*xbar is just the vector E, and that this vector is perpendicular
% to the columnspace of A.   

% Start with E.  We have already convinced ourselves of the fact that
% A' * E = 0. If you don't remember this, go back to the last statement 
% you executed.
% But E is just b - A*xbar
% Substituting for E, we get
% 		A' * (b - A*xbar) = 0
% rearranging, we get 
% 		A'*b = A'*A*xbar
% Remember, we are trying to solve for xbar. So we multiply both 
% sides of the equation by the inverse of (A'*A)
%		inv(A'*A) * A'*b 	= inv(A'*A) * (A'*A) * xbar
%		inv(A'*A) * A'*b 	= xbar

xbar = inv(A'*A) * A'*b 	% normal equation (see Strang, p. 156)

% The column vector describes the best fitting line to the data. 
% 1.2857 is the intercept and 0.5714 is the slope. 
 
figure(1)
fplot('1.2857 + 0.5714*x',[-3 3],'c')

% Now there is one thing left to appreciate before we extend this 
% to multiple regression. I promised that I would
% explain the projection matrix we used to make our pretty picture.
% It would be nice to be able to project a vector onto the columnspace of
% A.   Well, we already have one example:  A * xbar is the projection of
% b onto the columnspace.  Now look at what you get when you left-multiply the 
% equation for xbar by A.

% A * xbar = A * inv(A'*A) * A' * b

% This leads to a simple idea.  The seemingly messy string of matrices,
% A * inv(A'*A) * A', projects b onto the plane in our figure.

P = A * inv(A'*A) * A'

% Earlier in the tutorial, we used P to project the vertices of a square 
% onto the columnspace of A.  That's how we made the blue plane.

% So why is this so cool?  Because it gives us a very easy way to do
% multiple regression.  Suppose we have the ordered pairs,

data = [
     0     1
     1    53
     2    35
     3    -7
     4    87
     5   106
     6   250
     7   301
     8   346
     9   337
    10   256
    11   -25
    12  -421
]
		
figure(1)
hold off, plot(data(:,1),data(:,2),'o'),hold on

% We want to fit a 4th order polynomial.
% The A matrix has 5 columns

a = data(:,1);
b = data(:,2);
A = [ones(size(a)) a a.^2 a.^3 a.^4]

% We want the best solution for xbar in the equation A*xbar = b

xbar = inv(A'*A) * A'*b

% It's that simple.  We can get the fit back with
predvals = A*xbar
plot(a,predvals,'r*')

% and we can get a predicted curve by increasing the sampling density
a = linspace(min(a),max(a),100)';
A = [ones(size(a)) a a.^2 a.^3 a.^4];
predvals = A*xbar;
plot(a,predvals,'g:')


% If you want to go on, the next topic is weighted least squares
% Then, multiple regression with nondiagonal covariance matrix.
