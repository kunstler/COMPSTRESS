# Simple model of successional dynamics modified from Pacala and Rees (1998 Am. Nat.)

The model is based on two species, $E$ early successional species and $L$ late successional species. 
Three state for a cell $F$ free, and when occupied by a species, it can be in two states $._E$ early successional stage and $._L$ late successional stage. The dynamics of the system can be described by the following ode system of equations (which is a modified version of Pacala and Rees (1998), we will need to add $\alpha$ to model also a competition-colonisation trade-off).


\begin{equation}
\begin{array}{lcl}
\frac{dE_E}{dt} &=& c F (E_E + E_L) + c L_E (E_E + E_L) - \gamma E_E - \tau E_E , \\[6pt]
\frac{dE_L}{dt}  &=&   \gamma E_E - c E_E (L_E + L_L) - \tau E_L , \\[6pt]
\frac{dL_E}{dt}  &=&  c F (L_E + L_L) - c L_E (E_E + E_L) - \gamma L_E - \tau L_E , \\[6pt]
\frac{dL_L}{dt}  &=&  \gamma L_E + c E_E (L_E + L_L) - \tau L_L , \\[6pt]
F  &=&  1 - E_E - E_L - L_E - L_L ,
\end{array}
\end{equation}


with $\tau$ the mortality rate, $c$ the colonisation rate (constant across species), and $\gamma$ the succession rate.

Not so sure to find easily analytical solution (but numerical solving work fine ...)

At equilibrium $E_E^* + E_L^* + L_E^* + L_L^* = 1- \tau/c$ and $E_E^* + L_E^* = \tau/\gamma (E_L^* + L_L^* )$. Thus
\begin{equation}
E_E^* + L_E^* = \frac{\tau}{\tau + \gamma} (1 - \frac{\tau}{c}) = A
\end{equation}
and
\begin{equation}
E_L^* + L_L^* = \frac{\gamma}{\tau + \gamma} (1 - \frac{\tau}{c}) = 1 - A
\end{equation}

Then we need to work out $E_E^* + E_L^*$ and $L_E^* + L_L^*$.

Frome the first equation of the system sete to zero for equilibrium we get (considering that $F^* = \tau/c$):

\begin{equation}
0 = \tau  E_L^* + c L_E^* (E_E^* + E_L^*) - \gamma E_E^* , 
\end{equation}

Then as $L_E^* = A - E_E^*$ this gives:

\begin{equation}
0 = \tau  E_L^* + c (A - E_E^*) (E_E^* + E_L^*) - \gamma E_E^* , 
\end{equation}

Thus

\begin{equation}
E_L^* = E_E^* \frac{\gamma - c(A - E_E^*)}{\tau + c (A - E_E^*)}, 
\end{equation}

Then using the equation 3 of the system set to zero this gives:

\begin{equation}
0  =  \tau L_L^* - c L_E^* (E_E^* + E_L^*) - \gamma L_E^*  , 
\end{equation}

Thus 

\begin{equation}
L_L^*  =  \frac{1}{\tau} \bigg( c L_E^* (E_E^* + E_L^*) + \gamma L_E^* \bigg)  , 
\end{equation}

Replacing by the expression of $E_L^*$ and $L_E^*$ in function of $E_E^*$in gives:

\begin{equation}
L_L^*  =  \frac{(A - E_E^*) }{\tau} \Bigg( c E_E^*\big(\frac{ \tau +\gamma}{\tau + c (A - E_E^*)}\big) + \gamma \Bigg)  , 
\end{equation}

So we have

\begin{equation}
L_E^* = A - E_E^*
\end{equation}

\begin{equation}
E_L^* = E_E^* \frac{\gamma - c(A - E_E^*)}{\tau + c (A - E_E^*)}, 
\end{equation}

and

\begin{equation}
L_L^*  =  \frac{(A - E_E^*) }{\tau} \Bigg( c E_E^*\big(\frac{ \tau +\gamma}{\tau + c (A - E_E^*)}\big) + \gamma \Bigg)  , 
\end{equation}

In addition 

$E_E^* + E_L^* + L_E^* + L_L^* = 1- \tau/c$

We can thus derive an equation only based on $E_E^*$. 


$E_E^* + E_E^* \frac{\gamma - c(A - E_E^*)}{\tau + c (A - E_E^*)} + A - E_E^* + \frac{(A - E_E^*) }{\tau} \Bigg( c E_E^*\big(\frac{ \tau +\gamma}{\tau + c (A - E_E^*)}\big) + \gamma \Bigg) = 1- \tau/c$

And then express all quantity as fraction of $(\tau + c (A - E_E^*))$, assuming $\tau + c (A - E_E^*) >0$.

$\frac{1}{\tau + c (A - E_E^*)} \Bigg[E_E^* (\tau + c (A - E_E^*)) + E_E^* (\gamma - c(A - E_E^*)) +  (A - E_E^*)\big(\tau + c (A - E_E^*)\big) + \frac{(A - E_E^*) }{\tau} \Bigg( c E_E^* (\tau +\gamma) + \gamma (\tau + c (A - E_E^*)) \Bigg) \Bigg] = 1- \tau/c$

Thus


$E_E^* \bigg(\tau + c (A - E_E^*)\bigg) + E_E^* \bigg(\gamma - c(A -
E_E^*)\bigg) +  (A - E_E^*)\bigg(\tau + c (A - E_E^*)\bigg) +
\frac{(A - E_E^*) }{\tau} \Bigg( c E_E^* (\tau +\gamma) + \gamma
(\tau + c (A - E_E^*)) \Bigg)  = (1- \tau/c)\bigg(\tau + c (A -
E_E^*)\bigg)$


\begin{equation}
A(\tau + c A) + \frac{\gamma A}{\tau} (\tau + c A) - c A (1 + \frac{\gamma }{\tau}) E_E^* = \tau +
c A - \tau^2/c - \tau A +(\tau/c - 1) c E_E^*
\end{equation}



This gives a linear equation (strange no why does the degree 2 cancel out?):

$a_1 + a_2 E_E^* = 0$

with
\begin{equation}
a_1 = A (2 \tau+ \gamma - c) +c A^2 (1+\gamma / \tau) - \tau  +\tau^2/c
\end{equation}


\begin{equation}
a_2 = c (1 - \tau /c - A -\frac{A \gamma}{\tau})
\end{equation}


Thus
\begin{equation}
E_E^* = \frac{c \gamma \tau + \tau^3 -\tau^2 \gamma - \gamma^2 c}{c
\tau^2 + \gamma \tau - c \tau \gamma - c \gamma^2}
\end{equation}

For $E_E^* = 0$ $c \gamma \tau + \tau^3 -\tau^2 \gamma - \gamma^2 c =
0$
thus $\tau^3 + \gamma (\tau(c - \tau)) - c \gamma^2 = 0$

$\Delta = \tau^2 (c- \tau)^2 + 4 c \tau^3 = \tau^2 (c^2 +2c\tau
+\tau^2)= \tau^2 (c+\tau)^2$

$X_1 = \frac{-tau (c- \tau) - tau(c+\tau) }{-2c} = \tau$

There is a problem somewhere .......



## References

* Pacala, S.W. & Rees, M. (1998) Models suggesting field experiments to test two hypotheses explaining successional diversity. The American Naturalist, 152, 729–737.
