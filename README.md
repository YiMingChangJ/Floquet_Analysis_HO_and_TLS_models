# Floquet_Analysis_HO_and_TLS_models
\section{The Floquet states analysis}\label{sec:Floquet}
% \begin{equation}\label{eq: mechanical motion}
%     \mathbf{R}_{1}(t) = \alpha_{\rm m}\sin(\omega_{\rm m}t+\phi_m){\mathbf{x}}+x_0{\mathbf{x}}.
% \end{equation}
In this section, we discuss a 
method that allows one to
obtain 
analytical solutions to periodic linear differential equations,
specifically, using 
Floquet theory of quantum mechanical systems with a time-periodic Hamiltonian.
This will allow us to 
solve for the quasi-energy levels
in the presence of  a periodic dipole
term described from dipole spatial oscillation,
\begin{equation}\label{eq: mechanical motion}
    \mathbf{R}_{1}(t) = \alpha_{\rm m}\sin(\omega_{\rm m}t+\phi_m){\mathbf{x}}+x_0{\mathbf{x}}.
\end{equation}


% method to explain the new resonances of two coupled LOs with and without mechanical motion. We invoke the Floquet theory of quantum mechanical systems with a time-periodic Hamiltonian. 


Floquet engineering is a powerful tool in analysing periodic system
%or systematic calculations 
and also provides an intuitive understanding of time-driven system. One of the major strengths of Floquet theory, is that we can translate any periodic time-dependent Hamiltonian into an effective time-independent Hamiltonian represented by an infinity matrix~\cite{Shirley_1965}. Here we are interested in time-dependent Schr\"odinger equation as a differential equation in time
\begin{equation}\label{eq: TDSE}
    i\frac{\partial}{\partial t}\ket{\psi(t)} = H(t)\ket{\psi(t)},
\end{equation}
where the Hamiltonian inherits the time periodicity, $H(t) = H(t+nT)$ and $T = 2\pi/\omega_{\rm m}$, where $n$ is an integer, and $\omega_{\rm m}$ denotes the Floquet frequency or driving frequency. 

According to the Floquet theorem, there exists a complete set of solutions $\ket{\psi_\lambda(t)}$, indexed by $\lambda$~\cite{Shirley_1965,Breuer_2007}. It can be expressed in writing as 
\begin{equation}\label{eq: complete set}
    \ket{\psi_\lambda(t)} = e^{-i\varepsilon_\lambda t}\ket{\phi_\lambda (t)},
\end{equation}
where $\varepsilon_\lambda$ are real Floquet (time-independent) quasienergies or Floquet indices, and $\ket{\phi_\lambda(t)}$ are Floquet modes (states) with periodicity of $T$, $\ket{\phi_\lambda(t)}=\ket{\phi_\lambda(t+nT)}$. Substituting Eq.~\eqref{eq: complete set} into Eq.~\eqref{eq: TDSE}, 
%it immediately 
gives
\begin{equation}\label{eq: eigenvale equation}
    H_F\ket{\phi_\lambda(t)} = \varepsilon_\lambda\ket{\phi_\lambda(t)},
\end{equation}
where $H_F=H(t)-i\partial_t$ is the Floquet Hamiltonian operator, or so called the quasienergy operator. In addition, the completeness of the Floquet solution is expressed by the relation 
\begin{equation}
    \braket{\phi_\lambda|\phi_{\lambda'}} = \delta_{\lambda,\lambda'}
    ,\,\,\,\,\, 
    \sum_\lambda\ket{\phi_\lambda(t)}\bra{\phi_\lambda(t)}={I},
\end{equation}
which indicates every values of $t$ form a complete set of basis states. Thus, we can write any solution of time-dependent Schr\"odinger equation as a linear combination 
\begin{equation}\label{eq: general solution of TDSE}
    \ket{\psi(t)} = \sum_\lambda c_\lambda e^{-i\varepsilon_\lambda t}\ket{\phi_\lambda (t)},
\end{equation}
where $c_\lambda$ are time-independent coefficient.

\sh{Stuff below is two-level system specific, so should probably be done after introducing what a two level system is, right now you want to keep it general, before introducing the
dipole and dipole interaction representation of states for an atom}

For our system, we can write the state of the system as tensor product of two linear combinations of two-level system $\ket{\phi_\lambda(t)}=\ket{\alpha}\otimes\ket{\beta}$, where
\begin{equation}
    \begin{cases}
        \ket{\alpha} = c_{g_1,\lambda}(t)\ket{g_1} + c_{e_1,\lambda}(t)\ket{e_1},\\
        \ket{\beta} = c_{g_2,\lambda}(t)\ket{g_2} + c_{e_2,\lambda}(t)\ket{e_2}.
    \end{cases}
\end{equation}
 In addition, $c_{g,\lambda}(t)$ and $c_{e,\lambda}(t)$ are also periodic with period $T$, and $\ket{e}$ and $\ket{g}$ are excited and ground states. We then expand these time-dependent coefficients in Fourier series, $c_{g_1,\lambda}(t)=\sum_{m=-\infty}^\infty c_{g_1,\lambda}e^{im\omega_0 t}$, where $m$ is an integer, and insert into Eq.~\eqref{eq: eigenvale equation}:
\begin{equation}
    \sum_{m=-\infty}^\infty\sum_{\beta=e,g} (H_F-\varepsilon_\lambda) c_{\beta,\lambda}e^{im\omega_{\rm m} t}\ket{\phi_{\lambda'}} = 0. 
\end{equation}
Taking inner product with $\bra{\phi_\lambda}$, multiplying by $\exp{(-in\omega_0 T)}$, and integrating from 0 to $T$ with respect to time, we then obtain a matrix eigenvalue equation for the Floquet quasienergies:
\begin{multline}\label{eq: Floquet Hamiltonian matrix representation}
    \sum_{m,\beta}\left[\frac{1}{T}\int^T_0 dt H_{\lambda,\lambda'}e^{i(m-n)\omega_{\rm m} t} + n\omega_{\rm m}\delta_{\lambda,\lambda'}\delta_{m,n}\right]c^{(m)}_{\beta,\lambda}\\
    = \varepsilon_\lambda c^{(n)}_{\alpha,\lambda},
\end{multline}
where $H_{\lambda,\lambda'}$ = $\braket{\phi_\lambda|H(t)|\phi_\lambda'}$ and $\ket{\alpha},\ket{\beta} = \{\ket{e_j},\ket{g_j}\}$. The quantity in the square bracket of Eq.~\eqref{eq: Floquet Hamiltonian matrix representation} is the matrix representation of $H_F$, which is described by the row indexed by ($\lambda$, $n$) and column indexed by ($\lambda'$, $m$). This matrix can be constructed by restricting the number of integers $m$ and $n$, each combination of $(n, m)$ corresponds to a 4-by-4 block of elements match to the matrix components of a coupled two-level Hilbert space Hamiltonian, include $n\omega_0$ on the diagonal elements, thus we can calculate this matrix computationally. The index $\lambda$ labels eigenstates, and $m$ and $n$ are the Fourier mode indices. In order to determine the frequencies of the spectral lines present in the emitted spectrum, one can numerically calculate the eigenvalues of Eq.~\eqref{eq: Floquet Hamiltonian matrix representation} with respect to a specific truncation of integer values for $n$ and $m$, where $n, m=0,\pm,...,\pm l_{\rm max}$ and $l_{\rm max}$ is the order of the harmonics considered in the matrix. The total truncation is $N_l = 2l_{\rm max} +1$.%, which leads to $N_F = N_\lambda \cdot N_l$ and $N_\lambda$ represents the total states of two coupled TLSs. Therefore, the matrix representation $H_F$ is in size of $N_F\times N_F$


\subsection{Hamiltonian for Two Coupled Two-level System}

\sh{$e$ and $g$ also used above, so should be included here is introducing a coupled TLS example}

The atomic Hamiltonian describing a two-level system (TLS) is
\begin{equation}\label{eq: atomic Hamiltonian}
    H_{\rm atoms} = \frac{\hbar\omega_0}{2}{\sigma}_z = \frac{\hbar\omega_0}{2}(\ket{e}\bra{e}+\ket{g}\bra{g}\end{equation}
where $\hbar\omega_0$ is the energy difference between excited and ground states, ${\sigma}_z ={\sigma}^+{\sigma}^- - {\sigma}^-{\sigma}^+$ are Pauli matrices and $\omega_0$ is the transition frequency of the TLS. 

The exact {\it classical} equation that define the 
classical electric field generated by an idealized electric dipole located at the origin ($r=0$), is given by~\cite{jackson_classical_1999,novotny_hecht_2012} 
%\sh{add a reference, note this is also for free space, and d needs to be in the correct units - it is not the displacement, but a dipole moment (units of ?, so that $d \cdot E$ has units of energy. Explain how to  convert your e-nm to Debye (C-m) - e.g., 1 e-nm = $1.6 \times 10^{-19} \dot 10^{-9}$ c-m, so you need distance and electronic charge}
%sh{this is not really exact, but is the solution to an assumed harmonic field ${\bf E}(t)={\bf E}e^{-i\omega_0 t}$}
\begin{multline}\label{eq: classical E field}
\mathbf{E}(\mathbf{R}) = \frac{1}{4\pi\epsilon_0 \epsilon_b}\Bigg[k^2({\mathbf{n_R}}\times{\mathbf{d}})\times{\mathbf{n_R}}\frac{e^{ikR}}{R}+\\
\left[3({\mathbf{n_R}}\cdot{\mathbf{d}}){\mathbf{n_R}}-{\mathbf{d}}\right]\Bigg(\frac{1}{R^3}-\frac{ik}{R^2}\Bigg)e^{ikR}\Bigg],
\end{multline}
%\sh{You cannot use operators on the RHS and a classical field on the RHS, in this expression. In Eq. (9), ${\bf d}(t)$ is the classical time-depend dipole and $ {\bf r}$ is just a unit vector, so I woudl write this as ${\bf n}_{\bf r}$, then connect this classical expression to the interaction Hamiltonian needed in a dipole approximation. In these units, $d  = \epsilon_0 \epsilon \alpha E$, where  $\alpha$ is the polarizability volume (defined differently in different textbooks); you can also write equation (9) without the $\epsilon_0$ on teh bottom, and define teh dipole as $d  = \epsilon \alpha E$, but use the one consistent with the code}
where $\epsilon_b$ is
the electric permittivity
or dielectric constant
of 
the background material,
and we will consider free space, i.e., $\epsilon_b=n^2 =1$. 
Moreover, ${\mathbf{d}}={\alpha}\mathbf{E}$ 
% ${\mathbf{d}}_j = q \mathbf{r}_j$
is dipole moment, ${\bm \alpha}={\bf n_d} {\bf n_d}\alpha$ is the electric polarizability volume, $k = \omega_0/c$, 
%with $c$  the speed of the light,
and 
$R=|\mathbf{R}|$, and $\mathbf{n_r}=\mathbf{R}/R$.

For our system, we will consider two $s$-dipoles (i.e., $s$-polarized) that are coupled through the self-consistent electromagnetic field, so Dipole 1 interacts with the field generated by Dipole 2 (and vice versa). Thus,
as a first approximation,
and using a semiclassical approach, we can write our coupled dipoles interaction Hamiltonian with quantized dipoles as
\begin{align}\label{eq: interaction Hamiltonian}
    H_{\rm int} &= -{\mathbf{d}}_1\cdot{\mathbf{E}}_2(\mathbf{R}_{12}) =- \mathbf{d}_{1}{\sigma}_{x,1}\cdot\mathbf{E}_2(\mathbf{R}_{12}) \nonumber
    \\
      &=-\mathbf{d}_{1}\cdot\frac{1}{4\pi\epsilon_0}\Bigg[k_2^2({\mathbf{n_r}}\times\mathbf{d}_2)\times{\mathbf{n_r}}\frac{e^{ik_2R_{12}}}{R_{12}} \nonumber
    \\
    & \hspace{-0.8cm}
    +\left[3({\mathbf{n_r}}\cdot\mathbf{d}_2){\mathbf{n_r}}-\mathbf{d}_2\right]\Bigg(\frac{1}{R_{12}^3}-\frac{ik_2}{R_{12}^2}\Bigg)e^{ik_2R_{12}}\Bigg]
    {\sigma}_{x,1}\,{\sigma}_{x,2}\,, 
\end{align}
where ${\sigma}_{x,j}={\sigma}_{j}^++{\sigma}_{j}^-$, $\mathbf{R}_{12}=\mathbf{R}_1-\mathbf{R}_2$ is the displacement between centre of mass (COM) of Dipole 1 and Dipole 2, and $R_{12}=|\mathbf{R}_{12}|$. 

Recall, our mechanical motion formula Eq.~\eqref{eq: mechanical motion}, which is periodic time-dependent. Since the origin of the dipole moves as function of time, we can express the displacement between two dipoles as 
\begin{equation}\label{eq: displacement}
    \mathbf{R}_{12}(t) = \left[R_0 - \alpha_{\rm m}\sin(\omega_{\rm m} t)\right]\hat{x}.
\end{equation}
The approximation of the electric field in the near-field zone 
($|kR|\ll 1$, $e^{ikR}\approx 1$) is~\cite{novotny_hecht_2012}
\begin{equation}
    \mathbf{E}(\mathbf{R}_{12}(t))=\frac{1}{4\pi\varepsilon_0}[3({\mathbf{n_r}}\cdot\mathbf{d})\mathbf{n_r}-\mathbf{d}]\frac{1}{R_{12}^3(t)}.
\end{equation}
Thus, we can write the interaction Hamiltonian in the near field zone  as
% \begin{multline}\label{eq: t-depent interaction Hamiltonian}
% H_{\rm int}(t) = -\mathbf{d}_{1}\cdot\frac{1}{4\pi\varepsilon_0}\Bigg[k_2^2({\mathbf{r}}_{12}(t)\times\mathbf{d}_2)\times{\mathbf{r}}_{12}(t)\frac{e^{ik_2r_{12}(t)}}{r_{12}(t)}+
%     \\
% \left[3({\mathbf{r}}_{12}(t)\cdot\mathbf{d}_2){\mathbf{r}}_{12}(t)-\mathbf{d}_2\right]\Bigg(\frac{1}{r_{12}(t)^3}-\frac{ik_2}{r_{12}(t)^2}\Bigg)e^{ik_2r_{12}(t)}\Bigg]
% \\
% ({\sigma}_{x,1}{\sigma}_{x,2}).\\
% \end{multline}
% We are having $y$ polarized direction of dipoles and mechanical motion modify dipole horizontally ($x$-axis), thus $\mathbf{d}_j = y_j{y}$. After simplification, I have obtained
% \begin{multline}
%     H_{\rm int}(t) = -\frac{e^{ik_2r_{12}(t)}}{4\pi\varepsilon_0}y_1y_2({\sigma}_{x,1}{\sigma}_{x,2})\\
%     \left[k_2\frac{1}{r_{12}(t)}-\frac{1}{r_{12}^3(t)}+\frac{ik_2}{r_{12}^2(t)}\right]
% \end{multline}


\begin{equation}
    H_{\rm int, TLS}(t) = \frac{\mathbf{d}_1\mathbf{d}_2}{4\pi\epsilon_0}\frac{1}{R_{12}^3(t)} {\sigma}_{x,1}\, {\sigma}_{x,2}.
\end{equation}

Exploiting the geometric series
\begin{equation}
    \frac{1}{(1-r)^3} = 1 + 3r + 6r^2 + 10r^3 + 15r^4+...,
\end{equation}
with the term
\begin{equation}
    \frac{1}{R_{12}^3(t)} = \frac{1}{[R_0-{\alpha}_{\rm m} \sin(\omega_{\rm m} t)]^3} = \frac{1}{R_0^3\left[1-\tilde{\alpha}_{\rm m} \sin(\omega_{\rm m} t)
    \right]^3},
\end{equation}
 we can rewrite our interaction Hamiltonian as:
%
\begin{equation}
    H_{\rm int, TLS}(t) = \hbar g_{12}({\sigma}_{x,1}{\sigma}_{x,2})[1+3M+6M^2...],
\end{equation}
%
where $\hbar g_{12} = d_1 d_2/4\pi\epsilon_0 R_0^3$ is the static coupling strength, $M = \tilde{\alpha}_{\rm m} \sin(\omega_{\rm m} t)$ and $\tilde{\alpha}_{\rm m} = \alpha_{\rm m}/R_0$ is the normalized mechanical amplitude. Therefore, the full semi-classical time-dependent Hamiltonian of two coupled dipoles with mechanical motion can be written as
\begin{multline}\label{eq: Full Hamiltonian TLS}
    H_{\rm TLS}(t) \approx \frac{\hbar\omega_1}{2}{\sigma}_{z,1} +
    \frac{\hbar\omega_2}{2}{\sigma}_{z,2}+
    \hbar g_{12} {\sigma}_{x,1}{\sigma}_{x,2}
    \\
    +   \hbar g_{12} {\sigma}_{x,1}{\sigma}_{x,2} 
    \left[3\tilde{\alpha}_{\rm m}\sin(\omega_{\rm m}t)+6\tilde{\alpha}_{\rm m}^2\sin^2(\omega_{\rm m}t)\right.
    \\
    \left.+10\tilde{\alpha}_{\rm m}^3\sin^3(\omega_{\rm m}t) +15\tilde{\alpha}_{\rm m}^4\sin^4(\omega_{\rm m}t)\right],
\end{multline}
where we can separate the time-independent and time-dependent Hamiltonian as $H_{\rm TLS}(t)\rightarrow H_0' + H_{\rm AC}^F(t)$. The solution on this system returns four energy eigenvalues and states, ground, upper and lower polariton, and a two quanta excitation. 

\sh{Here I would just write:
$H= H_0 + H_{dd} + H^{\rm F}$
where $H_0 + H_{dd} = H_{\rm atoms}
+ \hbar g_{12} {\sigma}_{x,1}{\sigma}_{x,2}$ is time-independent.
}

% Are superradiant and subradiant states created by transition between resonance frequency ($2\omega_0$) and upper and lower polariton states?


\subsection{Fourier Approach with Dressed operator for Two Coupled Two-level Systems} \label{subsec: FA_TLS}

Next we 
%are able to 
transfer a time-dependent periodic problem to a time-independent one by inserting Eq.~\eqref{eq: Full Hamiltonian TLS} into Eq.~\eqref{eq: Floquet Hamiltonian matrix representation}, which can be solved by Fourier approach with a finite Floquet matrix. with a size of $N_F \times N_F$, where $N_F$ is equal to the multiplication number of truncation $N_l$ and number of energy states $N_\lambda$. The time-independent Hamiltonian can be diagonalized into a set of basis $\{E_\lambda,\ket{\phi_\lambda}\}$ ($H_0'\ket{\phi_\lambda} = E_\lambda\ket{\phi_\lambda}$). We can write the Hamiltonian matrix element as \red{we transfer the bare states (unperturbed emitter and field states) of the system to dressed states.}
\begin{equation}
    H_{\lambda\lambda'}(t)=E_\lambda\,\delta_{\lambda\lambda'}+A_{\lambda\lambda'}\left[ 3M+6M^2+...\right],
\end{equation}
where $\ket{\phi_\lambda} = \ket{\alpha}\otimes\ket{\beta}$ and $A_{\lambda\lambda'}= g_{12}\braket{\phi_\lambda|\sigma_{x,1}\sigma_{x,2}|\phi_{\lambda'}}$ with the matrix element of the Hermitian system operator. To ensure we have the correct states to mechanically perturb, we dress these two operators $\sigma_{x,1}$ and $\sigma_{x_2}$ with dressed states \sh{clarify - dress with dressed-sates sounds strange, dressed states are formed with dressing - e.g., the sub-radiant and super-radiant states are the dressed states of $H_0 + H_{\rm dd}$}, namely $\braket{\phi_\lambda|\sigma_{x,j}|\phi_{\lambda'}}$. Then, project our dressed-states operator in new basis to obtain new atom operator
\begin{equation}\label{eq: density matrix TLS}
    \tilde\sigma_j^+ = \sum_{i=0}^{N_\lambda} \sum_{k=i+1}^{N_\lambda} \braket{\phi_\lambda|\sigma_{x,j}|\phi_{\lambda'}}\ket{\phi_i} \bra{\phi_k},
\end{equation}
where $N_\lambda$ is number of states and $\tilde\sigma_j^-=(\tilde\sigma_j^+)^\dagger$. Thus, we can rewrite the Hamiltonian matrix element as
\begin{equation}\label{eq: matrix elements TLS}
    H_{\lambda\lambda'}(t)=E_\lambda\,\delta_{\lambda\lambda'}+g_{12}(\tilde\sigma_{x,1}\tilde\sigma_{x,2})\left[3M+6M^2+...\right].
\end{equation}

We then can substitute Eq.~\eqref{eq: matrix elements TLS} into the square bracket of Eq.~\eqref{eq: Floquet Hamiltonian matrix representation} to calculate Floquet Hamiltonian representation, by taking the integral of time-dependent elements with respect to time and switch sine to cosine 
%\ymc{(makes no difference, just having $\pm i$ in front of coefficients, but won't affect anything to the results)} 
to obtain our time-independent matrix,
\begin{equation}
\frac{1}{T}\int^T_0 e^{i(m-n)\omega_{\rm m} t} dt=
    \begin{cases}
        1, \text{   if   } m-n = 0,
        \\
        0,  \text{     if    } m \neq n.
    \end{cases}
\end{equation}

\begin{multline}
\frac{1}{T}\int^T_0 3\tilde{\alpha}_{\rm m}\cos(\omega_{\rm m} t)e^{i(m-n)\omega_{\rm m} t} dt=\\
    \begin{cases}
        \frac{3\tilde{\alpha}_{\rm m}}{2}, \text{   if   } m-n = \pm 1,
        \\
        0,  \text{     if    } m-n \neq \pm 1.
    \end{cases}
\end{multline}

\begin{multline}
\frac{1}{T}\int^T_0 6\tilde{\alpha}_{\rm m}^2\cos^2(\omega_{\rm m} t)e^{i(m-n)\omega_{\rm m} t} dt=\\
    \begin{cases}
        3\tilde{\alpha}_{\rm m}^2, \text{   if   } m-n = 0,
        \\
        \frac{3\tilde{\alpha}_{\rm m}^2}{2}, \text{   if   } m-n = \pm 2,
        \\
        0,  \text{     if    } m-n \neq 0,\,\pm 2.
    \end{cases}
\end{multline}

\begin{multline}
\frac{1}{T}\int^T_0 10\tilde{\alpha}_{\rm m}^3\cos^3(\omega_{\rm m} t)e^{i(m-n)\omega_{\rm m} t} dt=\\
    \begin{cases}
        \frac{15\tilde{\alpha}_{\rm m}^3}{4}, \text{   if   } m-n = \pm 1,
        \\
        \frac{5\tilde{\alpha}_{\rm m}^3}{4}, \text{   if   } m-n = \pm 3,
        \\
        0,  \text{     if    } m-n \neq \pm 1,\,\pm 3.
    \end{cases}
\end{multline}

\begin{multline}
\frac{1}{T}\int^T_0 15\tilde{\alpha}_{\rm m}^4\cos^4(\omega_{\rm m} t)e^{i(m-n)\omega_{\rm m} t} dt=\\
    \begin{cases}
        \frac{45\tilde{\alpha}_{\rm m}^4}{8}, \text{   if   } m-n = 0,
        \\
        \frac{15\tilde{\alpha}_{\rm m}^4}{4}, \text{   if   } m-n = \pm 2,
        \\
        \frac{15\tilde{\alpha}_{\rm m}^4}{16}, \text{   if   } m-n = \pm 4,
        \\
        0,  \text{     if    } m-n \neq 0,\,\pm 2,\,\pm 4.
    \end{cases}
\end{multline}
Therefore, the Floquet Hamiltonian representation matrix of TLS

\sh{you use TLS and TLA, I think using TLS is better and just say two level system}


model is given by
\begin{multline}\label{eq: Floquet_Hamiltonian_TLS}
    H_{F, \rm TLS} = \sum_{m,n}\sum_{\lambda,\lambda'}(E_\lambda+ n\omega_{\rm m})\delta_{\lambda,\lambda'}\delta_{m,n}
    \\
      +g_{12}(\tilde\sigma_{x,1}\tilde\sigma_{x,2})\left[\left(3\tilde{\alpha}_{\rm m}^2+\frac{45}{8}\tilde{\alpha}_{\rm m}^4\right)\delta_{m,n}\right.
     \\
     +\left(\frac{3\tilde{\alpha}_{\rm m}}{2}+\frac{15\tilde{\alpha}_{\rm m}^3}{4}\right)\delta_{m,n\pm 1}
     +\left(\frac{3\tilde{\alpha}_{\rm m}^2}{2}+\frac{15\tilde{\alpha}_{\rm m}^4}{4}\right)\delta_{m,n\pm 2}
     \\
     \left. +\frac{5\tilde{\alpha}_{\rm m}^3}{4}\delta_{m,n\pm 3} + \frac{15\tilde{\alpha}_{\rm m}^4}{16}\delta_{m,n\pm 4}\right].
\end{multline}
We have expanded the geometric series up to fourth order. 

\section{Analytical solution to two coupled classical oscillators and the quantum Hopfield model}

Recently, the
coupling between
a classical dipole and 
an oscillator representing a cavity mode was considered
without any RWA~\cite{hughes2023reconciling}
\sh{check RWA is defined}.
The solution for the exact polarizability, when coupled to the cavity mode, is
\begin{equation}
    \alpha(\omega) 
    = \frac{A_0\omega_0 }{\omega_0^2-\omega^2-
    (\omega_0/\omega_c)4 g^2\omega^2/(\omega_c^2-\omega^2)},
\end{equation}
where $4g^2 = d^2\omega_c/(2\epsilon_0V_{\rm eff} \epsilon_b)$,
$\omega_0$ is the resonance frequency of the TLS and $\omega_c$
is the resonance frequency of the cavity mode.
For on-resonance coupling, the pole solution for the resonances are
\begin{equation}\label{analytical_eigenfreq}
    \omega_\pm = \omega_0\sqrt{1+2\eta^2\pm 2\eta(1+\eta^2)^{1/2}},
\end{equation}
where $\eta = g/\omega_0$.

This solution was  shown to be identical to the poles of the quantum Hopfield model,
for the two excited resonances,
which in the dipole-gauge takes the form
($\hbar=1$)
\begin{equation}\label{eqn:Hop1}
    H_{\rm Hop} = \omega_c a^\dagger a + \omega_0 b^\dagger b + ig(a^{\dagger}-a)(b+b^\dagger)
    + D (b+b^\dagger)^2,
\end{equation}
where $\omega_c$ ($\omega_0$) is the cavity (atom) transition frequency,
and
%which includes a quadratic $P^2$ term,
%and we use 
$D=\eta g$ is the diamagnetic amplitude~\cite{Emary_2003,2002.02139}. 
For $\omega_0=\omega_c$, the eigenfrequency solutions are again 
$\omega_{\pm}$, and this is only the case when one includes the $D$ term. Note that the quadratic term here is known as $P^2$ term in a quantum 
dipole gauge, or as the $A^2$ term in the Coulomb gauge, which causes a spectral blueshift (and ensures unphysical behaviour when $\eta \geq 0.5$). 
Given the fact that 
when $g>0.1$, the quantum Hopfield model includes virtual photons in the ground state, and has squeezed states in vacuum (including the ground state), this is a somewhat remarkable result. 

Below, we will also show that 
these solutions coincide exactly
with the first two excited states of the coupled TLS
problem, when mechanical effects are off, when one defines
$\eta = g_{12}/\omega_0$
and considers
$\omega_1=\omega_2$.

\sh{Still need to think about explaining this, but it is an interesting fact, but I am puzzled a bit why there is only one Hopfield term if there are two TLSs, so it's not exactly the same quantum system - this model is for a quantum field dipole and an atom, so the physics is a bit different}


\section{Coupled Dipoles using a Quantized Harmonic Oscillators Model}
\subsection{Hamiltonian System of Quantized Harmonic Oscillators}
The atomic Hamiltonian of Quantum Harmonic Oscillators (QHO) describing two dipoles is
\begin{equation}\label{eq: atomic Hamiltonian HO}
    H_0 = \hbar\omega_0{a}^{\dagger}{a} + \hbar\omega_0 b^\dagger b,
\end{equation}
where ${a}^\dagger$, $b^\dagger$, and $b$, ${a}$ are creation and annihilation operators of Dipole 1 and 2. The derivation of QHO model is similar to TLS model, except $\sigma_x = a + a^\dagger$, and operators act on Fock states (number states of photons), which have infinity number of ladder states. Also, quantized dipole moment is $\mathbf{d}=q\hat{x} = q\sqrt{\frac{\hbar}{2m\omega_0}}({a}+{a}^\dagger)$. Thus, follow the same principles we can write the interaction Hamiltonian in the near field zone as
\begin{equation}\label{eq: H_int HO}
    H_{\rm int, QHO}(t) = \left(\frac{\hbar}{2m\omega_0}\right)\frac{q^2}{4\pi\epsilon_0}\frac{1}{R_{12}^3(t)}({\sigma}_{x,a}{\sigma}_{x,b}),
\end{equation}
where ${\sigma}_{x,a}=a+a^\dagger$ and ${\sigma}_{x,b} = b+b^\dagger$. In addition, by comparing the coefficients of interaction Hamiltonian of HO and TLS, it yields to $d_0^2\equiv q^2\hbar/2\omega_0m$, as expected. 
% \sh{where is $y_0$ defined here and if using a quantum dipole, you should use a hat, so that
% $\hat{\bf d}$, and how is that defined for the 
% TLS case? In both cases, the quantum dipole is simply
% ${\bf d} = {\bf d} \sigma_x$, and there is no need to define $d$ here unless you want to, as it will be described elsewhere, such as
% $q y_0$ for each dipole? The operators $a,b,\sigma$ are well known so we do not need to be hats there, but you can if you want}

Therefore, the full time-dependent Hamiltonian system of QHO up to fourth order of geometric series is given by 
\begin{multline}\label{eq: Full Hamiltonian HO}
    H_{\rm QHO} (t) \approx  \hbar\omega_0a^\dagger a+
    \hbar\omega_0 b^\dagger b+
    \hbar g_{12, \rm HO}({\sigma}_{x,a}{\sigma}_{x,b})
    \\
    [1+3\tilde{\alpha}_{\rm m}\sin(\omega_{\rm m}t)+6\tilde{\alpha}_{\rm m}^2\sin^2(\omega_{\rm m}t)
    \\
    +10\tilde{\alpha}_{\rm m}^3\sin^3(\omega_{\rm m}t) +15\tilde{\alpha}_{\rm m}^4\sin^4(\omega_{\rm m}t)],
\end{multline}
where $\hbar g_{12, HO} = q^2\hbar/6\pi m\omega_0\epsilon_0 R_0^3$. The solution on this system returns $N$ energy eigenvalues and states. 
\sh{Actually it returns an infinite number, and practical numerical simulations must check that the Hilbert space is large enough for both $a$ and $b$
to yield convergent results, for the states and energy regimes of interest. This is why the TLS operators are much better to use if the results are equivalent, with weak excitation, namely, under linear response. Bu tpractically we use 
$N= N_{\rm max}$, and define
$N_{\rm max}$ in terms of the number of states used for $a$ and $b$, soso
if $N_a$ and $N_b$, then
$N_{\rm max} = ?$}

\red{because photons are bosons, each state has an infinite ladder energy states that correspond to the number of photons that they contain.
}

\subsection{Fourier Approach with Dressed Operator for Harmonic Oscillator}
In order to avoid the unphysical results and recover the same linear upper and lower polariton states in the TLS Hamiltonian, we have added an $A^2$ term to our HO Hamiltonian system, namely the Hopfield model Eq.~\eqref{eqn:Hop1}, but with mechanical interactions,  so we
rewrite Eq.~\eqref{eq: Full Hamiltonian HO} as ($\hbar=1$)
\begin{multline}\label{eq: Full Hamiltonian HO, Hop1}
    H_{\rm QHO} (t) = \omega_0a^\dagger a+
    \omega_0 b^\dagger b +
    \\
    g_{12, \rm HO}({\sigma}_{x,a}{\sigma}_{x,b})[...]
    + (g_{12,\rm HO}^2/\omega_0){\sigma}_{x,b}^2[...]^2,
\end{multline}
where $[...] = \left[1+3\tilde{\alpha}_{\rm m}\sin(\omega_{\rm m}t)+6\tilde{\alpha}_{\rm m}^2\sin^2(\omega_{\rm m}t)\right.$
$\left.+10\tilde{\alpha}_{\rm m}^3\sin^3(\omega_{\rm m}t)+15\tilde{\alpha}_{\rm m}^4\sin^4(\omega_{\rm m}t)\right]$, is the time-dependent part.
Without time-dependence,
we have confirmed also numerically that the solution yields exactly the 
first two excited-state manifolds, are identical to the coupled TLS problem,
as we show below.
\sh{Confirm this is correct or clarify }


Following the same procedure as in Subsec.~\ref{subsec: FA_TLS}, we obtain the new dipole operator by projecting dressed-states operator in new basis, where we truncate the matrix operator from $N$-levels to a finite level ($N_{\rm trunc}$, a dressed density matrix). The new dipole operator is given by 
\begin{equation}\label{eq: density matrix HO}
    \tilde a^\dagger = \sum_{i=0}^{N_{\rm trunc}} \sum_{k=i+1}^{N_{\rm trunc}} \braket{j|(a + a^\dagger)|j'}\ket{\phi_i} \bra{\phi_k},
\end{equation}
where $\tilde a = (\tilde a^\dagger)^\dagger$. Therefore, the Floquet Hamiltonian representation matrix of HO model is given by
\begin{multline}\label{eq: Floquet_Hamiltonian_HO}
    H_{F, \rm HO} = \sum_{m,n}\sum_{j,j'}(E_j+ n\omega_{\rm m})\delta_{j,j'}\delta_{m,n}
    \\
      +g_{12, \rm HO}(\tilde\sigma_{x,a}\tilde\sigma_{x,b})\left[\left(3\tilde{\alpha}_{\rm m}^2+\frac{45}{8}\tilde{\alpha}_{\rm m}^4\right)\delta_{m,n}\right.
     \\
     +\left(\frac{3\tilde{\alpha}_{\rm m}}{2}+\frac{15\tilde{\alpha}_{\rm m}^3}{4}\right)\delta_{m,n\pm 1}
     +\left(\frac{3\tilde{\alpha}_{\rm m}^2}{2}+\frac{15\tilde{\alpha}_{\rm m}^4}{4}\right)\delta_{m,n\pm 2}
     \\
     \left. +\frac{5\tilde{\alpha}_{\rm m}^3}{4}\delta_{m,n\pm 3} + \frac{15\tilde{\alpha}_{\rm m}^4}{16}\delta_{m,n\pm 4}\right]
     \\
     +\frac{g_{12, \rm HO}^2}{\omega_0}\tilde\sigma_{x,b}^2\left[
     \left(\frac{21}{2}\tilde{\alpha}_{\rm m}^2+\frac{189}{4}\tilde{\alpha}_{\rm m}^4+\frac{175}{2}\tilde{\alpha}_{\rm m}^6+\frac{7875}{128}\tilde{\alpha}_{\rm m}^8\right)\delta_{m,n}\right.
     \\
     +\left(3\tilde{\alpha}_{\rm m}+21\tilde{\alpha}_{\rm m}^3+\frac{525}{8}\tilde{\alpha}_{\rm m}^5+\frac{2625}{32}\tilde{\alpha}_{\rm m}^7\right)\delta_{m,n\pm 1}
     \\
     +\left(\frac{21}{4}\tilde{\alpha}_{\rm m}^2+\frac{63}{2}\tilde{\alpha}_{\rm m}^4+\frac{525}{8}\tilde{\alpha}_{\rm m}^6+\frac{1575}{32}\tilde{\alpha}_{\rm m}^8\right)\delta_{m,n\pm 2}
     \\
     +\left(7\tilde{\alpha}_{\rm m}^3+\frac{525}{16}\tilde{\alpha}_{\rm m}^5+\frac{1575}{32}\tilde{\alpha}_{\rm m}^7\right)\delta_{m,n\pm 3}
     \\
     +\left(\frac{63}{8}\tilde{\alpha}_{\rm m}^4+\frac{105}{4}\tilde{\alpha}_{\rm m}^6+\frac{1575}{64}\tilde{\alpha}_{\rm m}^8\right)\delta_{m,n\pm 4}
     \\
     +\left(\frac{105}{16}\tilde{\alpha}_{\rm m}^5+\frac{525}{32}\tilde{\alpha}_{\rm m}^7\right)\delta_{m,n\pm 5}
     +\left(\frac{35}{8}\tilde{\alpha}_{\rm m}^6+\frac{225}{32}\tilde{\alpha}_{\rm m}^8\right)\delta_{m,n\pm 6}
     \\
     \left.
     +\left(\frac{75}{32}\tilde{\alpha}_{\rm m}^7\right)\delta_{m,n\pm 7}
     +\left(\frac{225}{256}\tilde{\alpha}_{\rm m}^8\right)\delta_{m,n\pm 8}
     \right],
\end{multline}
where $\tilde\sigma_{x,a} = \tilde a + \tilde a^\dagger$ and $\tilde\sigma_{x,b} = \tilde b + \tilde b^\dagger$, and the time-independent Hamiltonian is diagonalized into a set of basis $\{E_j,\ket{j}\}$. \red{since, $\left[1+3\tilde{\alpha}_{\rm m}\sin(\omega_{\rm m}t)+6\tilde{\alpha}_{\rm m}^2\sin^2(\omega_{\rm m}t)+...\right]^2$
not equal
$\left[3\tilde{\alpha}_{\rm m}\sin(\omega_{\rm m}t)+6\tilde{\alpha}_{\rm m}^2\sin^2(\omega_{\rm m}t)+...\right]^2$. I have to do the calculation of $[...]^2$ in hand, instead of numerically.}
