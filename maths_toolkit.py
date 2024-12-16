import numpy as np
from sympy import symbols, Eq, solve, Matrix, sin, cos, tan, cot, sec, csc, simplify

class MathsSolver:
    def __init__(self):
        pass

    @staticmethod
    def solve_linear_equations(equations, variables):
        """
        Solve a system of linear equations.
        :param equations: List of equations as strings.
        :param variables: List of variables as strings.
        :return: Solution as a dictionary.
        """
        print(f"Solving linear equations: {equations} with variables: {variables}")
        # Define variables as symbolic representations
        vars = symbols(variables)
        # Parse equations and convert to symbolic equality objects
        eqs = [Eq(eval(eq.split('=')[0], {str(var): var for var in vars}), eval(eq.split('=')[1], {str(var): var for var in vars})) for eq in equations]
        print(f"Parsed equations: {eqs}")
        # Solve the equations
        solution = solve(eqs, vars)
        print(f"Solution: {solution}")
        return solution

    @staticmethod
    def matrix_operations(matrix_a, matrix_b=None, operation="add"):
        """
        Perform matrix operations like addition, subtraction, multiplication.
        :param matrix_a: First matrix as a 2D list.
        :param matrix_b: Second matrix as a 2D list (optional).
        :param operation: Operation to perform - "add", "subtract", "multiply".
        :return: Resultant matrix as a 2D list.
        """
        print(f"Performing {operation} operation on matrices: {matrix_a} and {matrix_b}")
        # Convert input matrices to sympy Matrix objects
        mat_a = Matrix(matrix_a)
        mat_b = Matrix(matrix_b) if matrix_b else None

        # Perform the requested matrix operation
        if operation == "add":
            result = mat_a + mat_b
        elif operation == "subtract":
            result = mat_a - mat_b
        elif operation == "multiply":
            result = mat_a * mat_b
        else:
            raise ValueError("Invalid operation. Choose from 'add', 'subtract', 'multiply'.")
        print(f"Resultant matrix: {result}")
        return result.tolist()

    @staticmethod
    def find_locus(equation, x_var, y_var):
        """
        Find the locus of a given equation.
        :param equation: Equation as a string.
        :param x_var: Variable for x as a string.
        :param y_var: Variable for y as a string.
        :return: Rearranged equation representing the locus.
        """
        print(f"Finding locus for equation: {equation} with variables: {x_var}, {y_var}")
        # Define variables as symbolic representations
        x, y = symbols(f"{x_var} {y_var}")
        # Parse the equation and convert to a symbolic equality object
        eq = Eq(eval(equation.split('=')[0], {x_var: x, y_var: y}), eval(equation.split('=')[1], {x_var: x, y_var: y}))
        print(f"Parsed equation: {eq}")
        # Solve for the locus
        locus = solve(eq, y)
        print(f"Locus: {locus}")
        return locus

    @staticmethod
    def matrix_determinant(matrix):
        """
        Calculate the determinant of a matrix.
        :param matrix: Matrix as a 2D list.
        :return: Determinant as a float or integer.
        """
        print(f"Calculating determinant of matrix: {matrix}")
        # Convert input matrix to a sympy Matrix object
        mat = Matrix(matrix)
        # Compute the determinant
        determinant = mat.det()
        print(f"Determinant: {determinant}")
        return determinant

    @staticmethod
    def solve_quadratic(a, b, c):
        """
        Solve a quadratic equation ax^2 + bx + c = 0.
        :param a: Coefficient of x^2.
        :param b: Coefficient of x.
        :param c: Constant term.
        :return: Roots as a list.
        """
        print(f"Solving quadratic equation: {a}x^2 + {b}x + {c} = 0")
        # Define the variable as symbolic representation
        x = symbols('x')
        # Formulate the quadratic equation
        eq = Eq(a * x**2 + b * x + c, 0)
        print(f"Parsed equation: {eq}")
        # Solve the equation
        roots = solve(eq, x)
        print(f"Roots: {roots}")
        return roots

    @staticmethod
    def solve_trigonometric(equation, variable):
        """
        Solve trigonometric equations.
        :param equation: Equation as a string.
        :param variable: Variable as a string.
        :return: Solutions as a list.
        """
        print(f"Solving trigonometric equation: {equation} for variable: {variable}")
        # Define the variable as symbolic representation
        var = symbols(variable)
        # Parse the trigonometric equation and convert to a symbolic equality object
        eq = Eq(eval(equation.split('=')[0], {"sin": sin, "cos": cos, "tan": tan, "cot": cot, "sec": sec, "csc": csc, variable: var}),
                eval(equation.split('=')[1], {"sin": sin, "cos": cos, "tan": tan, "cot": cot, "sec": sec, "csc": csc, variable: var}))
        print(f"Parsed equation: {eq}")
        # Solve the trigonometric equation
        solutions = solve(eq, var)
        print(f"Solutions: {solutions}")
        return solutions

    @staticmethod
    def simplify_trigonometric(expression):
        """
        Simplify trigonometric expressions.
        :param expression: Expression as a string.
        :return: Simplified expression.
        """
        print(f"Simplifying trigonometric expression: {expression}")
        # Define the variable 'x' for use in trigonometric expressions
        x = symbols('x')
        # Parse the trigonometric expression
        expr = eval(expression, {"sin": sin, "cos": cos, "tan": tan, "cot": cot, "sec": sec, "csc": csc, "x": x})
        print(f"Parsed expression: {expr}")
        # Simplify the parsed expression
        simplified_expr = simplify(expr)
        print(f"Simplified expression: {simplified_expr}")
        return simplified_expr

# Example Usage
if __name__ == "__main__":
    solver = MathsSolver()

    # Example: Solve linear equations
    equations = ["2*x + y = 5", "x - y = 1"]
    variables = "x y"
    print("Linear Equations Solution:", solver.solve_linear_equations(equations, variables))

    # Example: Matrix operations
    mat_a = [[1, 2], [3, 4]]
    mat_b = [[5, 6], [7, 8]]
    print("Matrix Addition:", solver.matrix_operations(mat_a, mat_b, operation="add"))

    # Example: Find locus
    equation = "x**2 + y**2 - 25 = 0"
    print("Locus:", solver.find_locus(equation, "x", "y"))

    # Example: Determinant of a matrix
    matrix = [[2, 3], [1, 4]]
    print("Determinant:", solver.matrix_determinant(matrix))

    # Example: Solve quadratic equation
    print("Quadratic Roots:", solver.solve_quadratic(1, -3, 2))

    # Example: Solve trigonometric equation
    trig_eq = "sin(x) - 0.5 = 0"
    print("Trigonometric Solution:", solver.solve_trigonometric(trig_eq, "x"))

    # Example: Simplify trigonometric expression
    trig_expr = "sin(x)**2 + cos(x)**2"
    print("Simplified Expression:", solver.simplify_trigonometric(trig_expr))
