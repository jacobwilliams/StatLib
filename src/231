FUNCTION chisqnc(idf: nonnegint; nc, arg, eps: real; VAR density: real;
    VAR ifault: nonnegint): real;

{Algorithm AS 231  Appl. Statist. (1987) Vol.36, No.3}

{chisqnc evaluates the probability that a noncentral chi-squared 
    random variable is less than a given value, arg}

    CONST tol = -50;
          maxit = 500;
          lnrtpi2 = 0.22579135264473;

    VAR i, k: nonnegint;
        m: 0 .. maxit;
        ao, am, b1, eps2, hold, sum, dans, lans, pans, probability: real;
        converged: Boolean;

    PROCEDURE update(VAR j: nonnegint);
        BEGIN
        IF lans < tol THEN
            BEGIN
            lans := lans + ln(arg / j); dans := exp(lans)
            END
        ELSE dans := dans * arg / j;
        pans := pans - dans; j := j + 2
        END; {update}

    BEGIN
    IF (nc < 0.0) OR (arg < 0.0) OR (eps <= 0.0) THEN
        BEGIN
        chisqnc := -1.0; density := -1.0;
        ifault := 1
        END
    ELSE IF (arg = 0.0) AND (idf > 0) THEN
        BEGIN
        chisqnc := 0.0; density := 0.0;
        ifault := 0
        END
    ELSE
        BEGIN

        {Determine initial values}
        k := idf; b1 := 0.5 * nc;
        ao := exp(-b1); eps2 := eps / ao;
        IF odd(k) THEN
            BEGIN
            i := 1; lans := -0.5 * (arg + ln(arg)) - lnrtpi2;
            dans := exp(lans); pans := centnorm(sqrt(arg))
            END
        ELSE
            BEGIN
            i := 2; lans := -0.5 * arg;
            dans := exp(lans); pans := 1.0 - dans
            END;

        {Evaluate first term}
        IF k = 0 THEN
            BEGIN
            m := 1; k := 2;
            am := b1; sum := 1.0 / ao - 1.0 - am;
            density := am * dans; probability := 1.0 + am * pans
            END
        ELSE
            BEGIN
            m := 0; k := k - 1;
            am := 1.0; sum := 1.0 / ao - 1.0;
            WHILE i < k DO update(i);
            k := k + 1;
            density := dans; probability := pans
            END;

        {Evaluate successive terms of the expansion}
        REPEAT
            m := m + 1; am := b1 * am / m;
            update(k); sum := sum - am;
            density := density + am * dans; hold := am * pans;
            probability := probability + hold;
            converged := (pans * sum < eps2) AND (hold < eps2)
        UNTIL (m = maxit) OR converged;

        IF converged THEN ifault := 0 ELSE ifault := 2;
        density := 0.5 * ao * density;
        chisqnc := ao * probability
        END
    END; {of chisqnc}
