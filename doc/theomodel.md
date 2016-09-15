% Notes and ideas about theoretical model of tradeoff between competitive ability and stress tolerance
% Georges Kunstler, Irstea

# What are the key question that we can explore with this toy model?

Interplay between competition and abiotic stress is key to understand the changes in community structure along gradients and species ranges (Case et al. 2005).

A lot of studies have been focus on symmetric competition (see Case and Taper 2000, Norberg et al. 2012)

In forets communities competition is largely asymmetric BUT also driven by successional dynamics

Trade-off between competitive ability and climate stress tolerance

**Question:**

- can it produce realistic pattern of species distribution, diversity and species ranges?

# Simple model of successional dynamics modified from Pacala and Rees (1998 *Am. Nat.*)

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
E_E^* + L_E^* = \frac{\gamma}{\tau + \gamma} (1 - \frac{\tau}{c}) 
\end{equation}
and
\begin{equation}
E_L^* + L_L^* = \frac{\tau}{\tau + \gamma} (1 - \frac{\tau}{c}) 
\end{equation}

Then we need to work out $E_E^* + E_L^*$ and $L_E^* + L_L^*$.




## References

* Case, T.J. & Taper, M.L. (2000) Interspecific competition, environmental gradients, gene flow, and the coevolution of species’ borders. American Naturalist, 155, 583–605.
* Norberg, J., Urban, M.C., Vellend, M., Klausmeier, C.A. & Loeuille, N. (2012) Eco-evolutionary responses of biodiversity to climate change. Nature Climate Change, 2, 747–751.
* Pacala, S.W. & Rees, M. (1998) Models suggesting field experiments to test two hypotheses explaining successional diversity. The American Naturalist, 152, 729–737.
