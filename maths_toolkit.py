import math

class MathsToolkit:
    """A toolkit for intermediate-level Maths A and Maths B calculations."""

    # Algebra Functions
    @staticmethod
    def solve_quadratic(a, b, c):
        """Solve a quadratic equation ax^2 + bx + c = 0."""
        discriminant = b**2 - 4*a*c
        if discriminant < 0:
            return "No Real Roots"
        root1 = (-b + math.sqrt(discriminant)) / (2*a)
        root2 = (-b - math.sqrt(discriminant)) / (2*a)
        return root1, root2

    # Trigonometry Functions
    @staticmethod
    def sin_deg(angle):
        """Calculate sine of an angle in degrees."""
        return math.sin(math.radians(angle))

    @staticmethod
    def cos_deg(angle):
        """Calculate cosine of an angle in degrees."""
        return math.cos(math.radians(angle))

    # Calculus Functions
    @staticmethod
    def differentiate(poly_coeffs):
        """
        Differentiate a polynomial.
        poly_coeffs: List of coefficients [a_n, a_(n-1), ..., a_1, a_0].
        Returns: List of differentiated coefficients.
        """
        return [coeff * exp for exp, coeff in enumerate(poly_coeffs[::-1])][1:][::-1]

    @staticmethod
    def integrate(poly_coeffs):
        """
        Integrate a polynomial.
        poly_coeffs: List of coefficients [a_n, a_(n-1), ..., a_1, a_0].
        Returns: List of integrated coefficients.
        """
        return [coeff / (exp + 1) for exp, coeff in enumerate(poly_coeffs[::-1])][::-1] + [0]

    # Statistics Functions
    @staticmethod
    def mean(numbers):
        """Calculate the mean of a list of numbers."""
        return sum(numbers) / len(numbers) if numbers else 0

    @staticmethod
    def variance(numbers):
        """Calculate the variance of a list of numbers."""
        mean_val = MathsToolkit.mean(numbers)
        return sum((x - mean_val)**2 for x in numbers) / len(numbers) if numbers else 0

    @staticmethod
    def standard_deviation(numbers):
        """Calculate the standard deviation of a list of numbers."""
        return math.sqrt(MathsToolkit.variance(numbers))

    # Matrix Functions
    @staticmethod
    def matrix_add(matrix1, matrix2):
        """Add two matrices element-wise."""
        if len(matrix1) != len(matrix2) or any(len(row1) != len(row2) for row1, row2 in zip(matrix1, matrix2)):
            raise ValueError("Matrices must have the same dimensions for addition.")
        return [[val1 + val2 for val1, val2 in zip(row1, row2)] for row1, row2 in zip(matrix1, matrix2)]

    @staticmethod
    def matrix_multiply(matrix1, matrix2):
        """Multiply two matrices."""
        if len(matrix1[0]) != len(matrix2):
            raise ValueError("Number of columns in the first matrix must equal number of rows in the second matrix.")
        return [[sum(a * b for a, b in zip(row, col)) for col in zip(*matrix2)] for row in matrix1]

    @staticmethod
    def matrix_transpose(matrix):
        """Transpose a matrix."""
        return [list(row) for row in zip(*matrix)]
