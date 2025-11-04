def Newton(f, df, x0, epsilon, max_iter):
    iterations = [x0]  # Lista para armazenar todas as iterações
    if max_iter <= 0:
        raise ValueError("Número de iterações deve ser maior que zero.")
    
    print("k\tx\t\t\tf(x)")
    x = x0
    print("0\t%.8e\t%.8e" % (x, f(x)))
    
    for k in range(1, max_iter + 1):
        fx = f(x)
        if abs(fx) < epsilon:
            print("%d\t%.8e\t%.8e" % (k, x, fx))
            return x, iterations  # Retorna a raiz e as iterações
        
        dfx = df(x)
        if dfx == 0:
            raise ValueError("Derivada zero em x = %.8e" % x)
        
        x_new = x - fx / dfx
        iterations.append(x_new)  # Adiciona a nova iteração
        print("%d\t%.8e\t%.8e" % (k, x_new, f(x_new)))
        
        if abs(x_new - x) < epsilon:
            return x_new, iterations  # Retorna a raiz e as iterações
        
        x = x_new
    
    raise ValueError("Não convergiu após %d iterações (último x: %.8e)" % (max_iter, x))