## for the varE calc, we need <H^2>. This computation can be simplified using the anticommutator, which has been implemented similarly to the commutator of the PP library
## the varE calc is not the bottleneck of the code (rather the VQE itself).

using PauliPropagation

"""
    anticommutator(pstr1::PauliStringType, pstr2::PauliStringType)

Calculate the anticommutator of two integer Pauli strings.
Returns (new_pstr, coeff). The coefficient is zero if the Pauli strings anticommute.
"""
function anticommutator(pstr1::PauliStringType, pstr2::PauliStringType)
    if commutes(pstr1, pstr2)
        # For commuting pairs, pauliprod returns a real sign in {+1,-1}
        new_pstr, total_sign = pauliprod(pstr1, pstr2)  # total_sign ∈ {+1,-1} (ComplexF64 but real)
        return new_pstr, 2.0 * total_sign
    else
        # {A,B}=0 if they anticommute
        return identitylike(pstr1), ComplexF64(0.0)
    end
end


# copy-paste of the commutator wrapper of PP
"""
    anticommutator(pstr1::PauliString, pstr2::PauliString)
"""
function anticommutator(pstr1::PauliString, pstr2::PauliString)
    new_pstr, new_coeff = anticommutator(pstr1.term, pstr2.term)
    return PauliString(pstr1.nqubits, new_pstr, new_coeff)
end

"""
    anticommutator(psum::PauliSum, pstr::PauliString)
    anticommutator(pstr::PauliString, psum::PauliSum)
"""
anticommutator(psum::PauliSum, pstr::PauliString) = anticommutator(psum, PauliSum(pstr))
anticommutator(pstr::PauliString, psum::PauliSum) = anticommutator(PauliSum(pstr), psum)

"""
    anticommutator(psum1::PauliSum, psum2::PauliSum)
"""
function anticommutator(psum1::PauliSum, psum2::PauliSum)
    new_pstr_dict = anticommutator(psum1.terms, psum2.terms)
    return PauliSum(psum1.nqubits, new_pstr_dict)
end


# identical structure to the commutator Dict method, but flipped the outer `if`
function anticommutator(psum1::Dict{TT,CT1}, psum2::Dict{TT,CT2}) where {TT,CT1,CT2}
    new_pauli_dict = Dict{TT,ComplexF64}()

    for (pauli1, coeff1) in psum1, (pauli2, coeff2) in psum2
        if commutes(pauli1, pauli2)               
            new_pstr, s = anticommutator(pauli1, pauli2)  # s ∈ {±2, 0}
            if !iszero(s)                           # skip identity bookkeeping when zero
                new_pauli_dict[new_pstr] =
                    get(new_pauli_dict, new_pstr, 0.0 + 0.0im) + s * coeff1 * coeff2
            end
        end
    end

    # delete the pauli strings with zero coeffs
    for k in collect(keys(new_pauli_dict))
        if iszero(new_pauli_dict[k]) || abs(new_pauli_dict[k]) ≈ 0.0
            delete!(new_pauli_dict, k)
        end
    end

    return new_pauli_dict
end