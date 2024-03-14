
:- use_module(library(apply)).
:- use_module(library(lists)).
:- use_module(library(aggregate)).

%%% Defining diseases %%%
    
1/3::disease("x"); 1/3::disease("y"); 1/3::disease("auto").

%%% X-linked alleles
allele("xx", "aa").
allele("xx", "Aa").
allele("xx", "AA").

allele("xy", "a-").
allele("xy", "A-").

%%% Y-linked alleles
allele("-", "-").
allele("y", "a").
allele("y", "A").

%%% Autosomal alleles
allele("auto", "aa").
allele("auto", "Aa").
allele("auto", "AA").

%%% Defining alleles %%%

f_allele("xx", X) :- allele("xx", X).
f_allele("-", X) :- allele("-", X).
f_allele("auto", X) :- allele("auto", X).
    
m_allele("xy", X) :- allele("xy", X).
m_allele("y", X) :- allele("y", X).    
m_allele("auto", X) :- allele("auto", X).

%%% Possible chromosomes depending on diseases %%%

allele_disease("xy", "x").
allele_disease("xx", "x").
allele_disease("auto", "auto").
allele_disease("-", "y").    
allele_disease("y", "y").   

%% possible_allele(A) :- allele(X, A), allele_disease(X, D), disease(D).

%%% Defining dominance / recessiveness / penetrance

show("A", "AA"). % by default dominant allele
show("a", "aa"). % by default recessive allele

relative_penetrance("A", 1). % put 1/2 wherein co-dominance
relative_penetrance("a", 0). % put 1/2 wherein co-dominance

show("A",  "Aa") :- relative_penetrance("A", R1), relative_penetrance("a", R2), R1  > R2.
show("Aa", "Aa") :- relative_penetrance("A", R1), relative_penetrance("a", R2), R1 =:= R2. 
show("a",  "Aa") :- relative_penetrance("A", R1), relative_penetrance("a", R2), R1  < R2.
    
show("A", "A-").
show("a", "a-").
show("-", "-").
show("A", "A").
show("a", "a").

penetrance("-", 1).
penetrance("a", 1).
penetrance("A", 1).
penetrance("Aa", 1).  

%%%%% X-linked alleles initial probabilities %%%%%

P::m_allele_hw(I, "xy", "a-"); (1-P)::m_allele_hw(I, "xy", "A-") :- 
                    disease("x"), prevalence("a", P), generation(I).
                    
P**2::f_allele_hw(I, "xx", "aa"); 
    (2*P-2*P**2)::f_allele_hw(I, "xx", "Aa"); 
        (1-P)**2::f_allele_hw(I, "xx", "AA") :- 
                    disease("x"), prevalence("a", P), generation(I).
                    
%%%%% Y-linked alleles initial probabilities %%%%%

P::m_allele_hw(I, "y", "a"); (1-P)::m_allele_hw(I, "y", "A") :- 
                    disease("y"), prevalence("a", P), generation(I).
                    
f_allele_hw(I, "-", "-") :- disease("y"), generation(I).
    
%%%%% Autosomal alleles initial probabilities %%%%%
                    
P**2::m_allele_hw(I, "auto", "aa"); 
        (2*P-2*P**2)::m_allele_hw(I, "auto", "Aa"); 
            (1-P)**2::m_allele_hw(I, "auto", "AA") :- 
                    disease("auto"), prevalence("a", P), generation(I).
                            
P**2::f_allele_hw(I, "auto", "aa"); 
        (2*P-2*P**2)::f_allele_hw(I, "auto", "Aa"); 
            (1-P)**2::f_allele_hw(I, "auto", "AA") :- 
                    disease("auto"), prevalence("a", P), generation(I). 
    
%%%% Predicate decides gender of family descendant at generation I %%%%

P::m_family_descendant(I); (1-P)::f_family_descendant(I) :- generation(I), I > 0, P is 1/2.
                    
%%%% Apply hw principle to generation zero & onwards %%%%
                            
m_allele(0, C, A) :- m_allele_hw(0, C, A), allele(C, A), allele_disease(C, D), disease(D).
f_allele(0, C, A) :- f_allele_hw(0, C, A), allele(C, A), allele_disease(C, D), disease(D).
                                     
%%% Calculating offspring genotype wherein autosomal disorder

1/2::offspring_genotype(II, "auto", "Aa") ; 
1/4::offspring_genotype(II, "auto", "aa") ; 
1/4::offspring_genotype(II, "auto", "AA") :-  
            m_allele(I, "auto", "Aa"), 
            f_allele(I, "auto", "Aa"),
            disease("auto"), 
            generation(I), I >= 0, II is I+1.
        
1/2::offspring_genotype(II, "auto", "Aa") ; 
1/2::offspring_genotype(II, "auto", "aa") :- 
            m_allele(I, "auto", "Aa"), 
            f_allele(I, "auto", "aa"), 
            disease("auto"),
            generation(I), I >= 0, II is I+1.
            
1/2::offspring_genotype(II, "auto", "Aa") ;
1/2::offspring_genotype(II, "auto", "AA") :- 
            m_allele(I, "auto", "Aa"), 
            f_allele(I, "auto", "AA"),
            disease("auto"), 
            generation(I), I >= 0, II is I+1.
            
1/2::offspring_genotype(II, "auto", "Aa") ;
1/2::offspring_genotype(II, "auto", "aa") :- 
            m_allele(I, "auto", "aa"), 
            f_allele(I, "auto", "Aa"), 
            disease("auto"),
            generation(I), I >= 0, II is I+1.
            
1/2::offspring_genotype(II, "auto", "Aa") ; 
1/2::offspring_genotype(II, "auto", "AA") :-
            m_allele(I, "auto", "AA"), 
            f_allele(I, "auto", "Aa"), 
            disease("auto"),
            generation(I), I >= 0, II is I+1.

offspring_genotype(II, "auto", "aa") :- m_allele(I, "auto", "aa"), f_allele(I, "auto", "aa"), disease("auto"), 
    generation(I), I >= 0, II is I+1.
offspring_genotype(II, "auto", "AA") :- m_allele(I, "auto", "AA"), f_allele(I, "auto", "AA"), disease("auto"),
    generation(I), I >= 0, II is I+1.
            
offspring_genotype(II, "auto", "Aa") :- m_allele(I, "auto", "aa"), f_allele(I, "auto", "AA"), disease("auto"),
    generation(I), I >= 0, II is I+1.
offspring_genotype(II, "auto", "Aa") :- m_allele(I, "auto", "AA"), f_allele(I, "auto", "aa"), disease("auto"),
    generation(I), I >= 0, II is I+1.
    
%% inheritance respects hw, independant of sex-linked heritance since on autosome
    
m_allele(I, "auto", A) :- f_family_descendant(I), m_allele_hw(I, "auto", A), allele("auto", A), disease("auto"), generation(I), I > 0.  
f_allele(I, "auto", A) :- m_family_descendant(I), f_allele_hw(I, "auto", A), allele("auto", A), disease("auto"), generation(I), I > 0.

m_allele(I, "auto", A) :- m_inherit_genotype(I, "auto", A), allele("auto", A), disease("auto"), generation(I), I > 0.  
f_allele(I, "auto", A) :- f_inherit_genotype(I, "auto", A), allele("auto", A), disease("auto"), generation(I), I > 0.

%% inheritance independant of sex-linked heritance since on autosome

inherit_genotype(I, "auto", A) :- offspring_genotype(I, "auto", A), disease("auto"), generation(I), I > 0.
    
m_inherit_genotype(I, "auto", A) :- inherit_genotype(I, "auto", A), m_family_descendant(I).
f_inherit_genotype(I, "auto", A) :- inherit_genotype(I, "auto", A), f_family_descendant(I).
    
genotype_new(I, "auto", A) :- m_allele(I, "auto", A), f_inherit_genotype(I, "auto", A), f_family_descendant(I).
genotype_new(I, "auto", A) :- f_allele(I, "auto", A), m_inherit_genotype(I, "auto", A), m_family_descendant(I).

%%% Calculating offspring genotype wherein x-linked disorder

1/2::offspring_genotype(II, "xx", "AA") ; 
1/2::offspring_genotype(II, "xx", "Aa") :- 
        m_allele(I, "xy", "A-"), 
        f_allele(I, "xx", "Aa"),
        f_family_descendant(II),
        disease("x"),
        generation(I), I >= 0, II is I+1.
        
1/2::offspring_genotype(II, "xy", "A-") ;
1/2::offspring_genotype(II, "xy", "a-") :- 
        m_allele(I, "xy", "A-"), 
        f_allele(I, "xx", "Aa"),
        m_family_descendant(II),
        disease("x"),
        generation(I), I >= 0, II is I+1.        

1/2::offspring_genotype(II, "xx", "Aa") ;
1/2::offspring_genotype(II, "xx", "aa") :- 
            m_allele(I, "xy", "a-"), 
            f_allele(I, "xx", "Aa"),
            f_family_descendant(II),
            disease("x"),
            generation(I), I >= 0, II is I+1.        
        
1/2::offspring_genotype(II, "xy", "A-") ;
1/2::offspring_genotype(II, "xy", "a-") :- 
            m_allele(I, "xy", "a-"), 
            f_allele(I, "xx", "Aa"), 
            m_family_descendant(II),
            disease("x"),
            generation(I), I >= 0, II is I+1.

offspring_genotype(II, "xx", "AA") :- m_allele(I, "xy", "A-"), f_allele(I, "xx", "AA"), 
            f_family_descendant(II), disease("x"), generation(I), I >= 0, II is I+1.
        
offspring_genotype(II, "xy", "A-") :- m_allele(I, "xy", "A-"), f_allele(I, "xx", "AA"), 
            m_family_descendant(II), disease("x"), generation(I), I >= 0, II is I+1.

offspring_genotype(II, "xx", "aa") :- m_allele(I, "xy", "a-"), f_allele(I, "xx", "aa"), 
            f_family_descendant(II), disease("x"), generation(I), I >= 0, II is I+1.

offspring_genotype(II, "xy", "a-") :- m_allele(I, "xy", "a-"), f_allele(I, "xx", "aa"), 
            m_family_descendant(II), disease("x"), generation(I), I >= 0, II is I+1.

offspring_genotype(II, "xx", "Aa") :- m_allele(I, "xy", "A-"), f_allele(I, "xx", "aa"), 
            f_family_descendant(II), disease("x"), generation(I), I >= 0, II is I+1.
            
offspring_genotype(II, "xy", "a-") :- m_allele(I, "xy", "A-"), f_allele(I, "xx", "aa"), 
            m_family_descendant(II), disease("x"), generation(I), I >= 0, II is I+1.

offspring_genotype(II, "xx", "Aa") :- m_allele(I, "xy", "a-"), f_allele(I, "xx", "AA"), 
            f_family_descendant(II), disease("x"), generation(I), I >= 0, II is I+1.
            
offspring_genotype(II, "xy", "A-") :- m_allele(I, "xy", "a-"), f_allele(I, "xx", "AA"), 
            m_family_descendant(II), disease("x"), generation(I), I >= 0, II is I+1.
        
% Hardy–Weinberg principle
m_allele(I, "xy", A) :- f_family_descendant(I), m_allele_hw(I, "xy", A), allele("xy", A), disease("x"), generation(I), I > 0.  
f_allele(I, "xx", A) :- m_family_descendant(I), f_allele_hw(I, "xx", A), allele("xx", A), disease("x"), generation(I), I > 0.
    
m_allele(I, "xy", A) :- m_inherit_genotype(I, "xy", A), allele("xy", A), disease("x"), generation(I), I > 0.  
f_allele(I, "xx", A) :- f_inherit_genotype(I, "xx", A), allele("xx", A), disease("x"), generation(I), I > 0.   
    
% Probability inheritance wth equal M/F ratio
m_inherit_genotype(I, "xy", A) :- offspring_genotype(I, "xy", A), m_family_descendant(I).
f_inherit_genotype(I, "xx", A) :- offspring_genotype(I, "xx", A), f_family_descendant(I).
    
genotype_new(I, "xy", A) :- m_allele(I, "xy", A), f_inherit_genotype(I, "xy", A), f_family_descendant(I).
genotype_new(I, "xx", A) :- f_allele(I, "xx", A), m_inherit_genotype(I, "xx", A), m_family_descendant(I).

%%% Calculating offspring genotype wherein y-linked disorder        

offspring_genotype(II, "y", "a") :- 
            m_allele(I, "y", "a"),
            m_family_descendant(II),
            disease("y"),  
            generation(I), I >= 0, II is I+1. 
            
offspring_genotype(II, "y", "A") :- 
            m_allele(I, "y", "A"),
            m_family_descendant(II),  
            disease("y"),
            generation(I), I >= 0, II is I+1. 
    
offspring_genotype(II, "-", "-") :- 
            f_allele(I, "-", "-"),
            f_family_descendant(II),
            disease("y"),
            generation(I), I >= 0, II is I+1.         
    
m_allele(I, "y", A) :- f_family_descendant(I), m_allele_hw(I, "y", A), allele("y", A), disease("y"), generation(I), I > 0. 
f_allele(I, "-", A) :- m_family_descendant(I), f_allele_hw(I, "-", A), allele("-", A), disease("y"), generation(I), I > 0.
    
m_allele(I, "y", A) :- m_inherit_genotype(I, "y", A), allele("y", A), disease("y"), generation(I), I > 0.  
f_allele(I, "-", A) :- f_inherit_genotype(I, "-", A), allele("-", A), disease("y"), generation(I), I > 0.       
    
m_inherit_genotype(I, "y", A) :- offspring_genotype(I, "y", A), m_family_descendant(I).
f_inherit_genotype(I, "-", A) :- offspring_genotype(I, "-", A), f_family_descendant(I).
    
genotype_new(I, "y", A) :- m_allele(I, "y", A), f_inherit_genotype(I, "y", A), f_family_descendant(I).
genotype_new(I, "-", A) :- f_allele(I, "-", A), m_inherit_genotype(I, "-", A), m_family_descendant(I).
                       
%%% Carrying allele meta-relationship                        
                        
m_carry(I, A) :- m_allele(I, C, A), allele(C, A), allele_disease(C, D), disease(D), generation(I), I>=0.
f_carry(I, A) :- f_allele(I, C, A), allele(C, A), allele_disease(C, D), disease(D), generation(I), I>=0.
    
%%% Showing phenotype relationship

P::m_show(I, Ph) :- show(Ph, A), m_allele(I, C, A), allele_disease(C, D), disease(D), penetrance(Ph, P), generation(I), I>=0.
P::f_show(I, Ph) :- show(Ph, A), f_allele(I, C, A), allele_disease(C, D), disease(D), penetrance(Ph, P), generation(I), I>=0.

