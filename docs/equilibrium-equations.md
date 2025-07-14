# Equilibrium Equations

## General gas phase reaction:

$$aA + bB \leftrightarrow rR + sS$$

- Find reaction conversion under any $T$ and $P$ conditions

1. First calculate equilibrium constants $K$  
2. Use $K$ to calculate position of equilibrium (and therefore reaction extent $\chi$)

---

## Equilibrium Constant

Equilibrium constant equation:

$$K = \exp \left( \frac{-\Delta G}{RT} \right)$$

- $\Delta G$ depends on reaction conditions

$$\Delta G = \Delta H - T \Delta S$$

- $\Delta H$ and $\Delta S$ depend on temperature; use Hess cycle to find values for any $T$:

Notes:  
- Index $i$ refers to reactant species $i$ in reaction system  
- Index $j$ refers to product species $j$ in reaction system  
- $\gamma$ is the stoichiometric coefficient of a given component in the reaction system

Temperature adjusted enthalpy change:

$$\Delta H_T = q_1 + \Delta H^{\ominus} + q_2$$

$$\Delta H_T = \sum_{i} \int_{T}^{298} \gamma_i C_{p_i}(T) \, dT + \Delta H^{\ominus} + \sum_{j} \int_{298}^{T} \gamma_j C_{p_j}(T) \, dT$$

Temperature adjusted entropy change:

$$\Delta S_T = q_1 + \Delta S^{\ominus} + q_2$$

$$\Delta S_T = \sum_{i} \int_{T}^{298} \gamma_i \frac{C_{p_i}(T)}{T} \, dT + \Delta S^{\ominus} + \sum_{j} \int_{298}^{T} \gamma_j \frac{C_{p_j}(T)}{T} \, dT$$

Temperature adjusted free energy change:

$$\Delta G_T = \Delta H_T - T \Delta S_T$$

$$K = \exp \left( \frac{-\Delta G_T}{RT} \right)$$

- $K$ is dimensionless  
- $K$ valid for temperature $T$ at $1 \text{atm}$ due to standard ($\ominus$) enthalpy and entropy used to calculate adjusted values  
- Need to adjust for pressure $P$

---

## Reaction Quotient

Notes:  
- Reaction extent will be characterised by dimensionless quantity $\chi$  
- $I$ is the initial number of moles of a given component in the reaction system at time $t=0$  
- Aim is to solve for $\chi$

Consider gas phase reaction $aA + bB \leftrightarrow rR + sS$:

Reaction material balance:

| Species      | $aA$            | $bB$            | $rR$            | $sS$            |
|--------------|-------------------|-------------------|-------------------|-------------------|
| Initial moles| $I_A$           | $I_B$           | $I_R$           | $I_S$           |
| Change       | $-a \chi$        | $-b \chi$       | $+r \chi$       | $+s \chi$       |
| Final moles  | $I_A - a \chi$   | $I_B - b \chi$  | $I_R + r \chi$  | $I_S + s \chi$  |

Reaction quotient expression for ideal gas phase reaction:

$$Q = \frac{P_R^r P_S^s}{P_A^a P_B^b} = \frac{y_R^r \cdot y_S^s}{y_A^a \cdot y_B^b} \cdot P^{(r + s) - (a + b)}$$

- $y$ is gaseous mole fraction of a component

Mole fractions are related to total moles:

$$y_A = \frac{n_A}{n_{\text{total}}}$$

Where $n$ is the number of moles of a component at equilibrium (time $t = \infty$).

The terms $n_x$ and $n_{\text{total}}$ contain $\chi$ so all $n$ terms must be written in terms of $\chi$:

$$n_{\text{total}} = \sum_i (I_i - \gamma_i \chi) + \sum_j (I_j + \gamma_j \chi)$$

$$n_x = I_x \pm \gamma_x \chi$$

Rewrite reaction quotient $Q$ in terms of $\chi$:

$$Q = \frac{(I_R + \gamma_R \chi)^r \cdot (I_S + \gamma_S \chi)^s}{(I_A - \gamma_A \chi)^a \cdot (I_B - \gamma_B \chi)^b} \cdot \left(\frac{1}{\sum_i (I_i - \gamma_i \chi) + \sum_j (I_j + \gamma_j \chi)}\right)^{\Delta \gamma} \cdot P^{\Delta \gamma}$$

- $\Delta \gamma = (r + s) - (a + b)$

Final step is to equate $Q$ to $K$ at equilibrium:

$$Q - K = 0$$

Solve this equation for $\chi$ to find reaction extent. Only valid if $Q$ is also dimensionless.

### Pressure dependance

Define $P$ as pressure ratio of reaction pressure to standard (atmospheric) pressure:
$$
P = \frac{P_{\text{reac}}}{P_{\text{atm}}}
$$

So $P^{\Delta \gamma}$ is dimensionless, allowing $Q$ to be dimensionless.

---

## Generalised Form

Reaction quotient:

$$Q = \frac{\prod_j (I_j + \gamma_j \chi)^{\gamma_j}}{\prod_i (I_i - \gamma_i \chi)^{\gamma_i}} \cdot \left(\frac{1}{\sum_i (I_i - \gamma_i \chi) + \sum_j (I_j + \gamma_j \chi)}\right)^{\Delta \gamma} \cdot P^{\Delta \gamma}$$
