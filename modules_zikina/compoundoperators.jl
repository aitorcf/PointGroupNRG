include( "symbols.jl" )
include( "numericals.jl" )

# ############# #
# MISCELLANEOUS #
# ............. #
function show_operator_action( o::Operator )
    for s in o.basis.states 
        println( string(s) , " ==> " , string(o*State(s,o.basis)) )
    end
end

function show_operator_average( o::Operator )
    for s in o.basis.states 
        println( string(s) , " ==> " , string(State(s,o.basis)*o*State(s,o.basis)) )
    end
end


# #### #
# SPIN #
# .... #
function pauli_z( snospin::Tuple , hilbert::HilbertSpace , basis="canonical" )
    basis=="canonical" && (basis = CanonicalBasis( hilbert ))
    op = Operator( SymbolCreationOperator((snospin...,"u")) , basis )*
         Operator( SymbolAnnihilationOperator((snospin...,"u")) , basis )-
         Operator( SymbolCreationOperator((snospin...,"d")) , basis )*
         Operator( SymbolAnnihilationOperator((snospin...,"d")) , basis )
    return op
end
spin_z( snospin::Tuple , hilbert::HilbertSpace , basis="canonical" ) =
    0.5 * pauli_z( snospin , hilbert )

function pauli_x( snospin::Tuple , hilbert::HilbertSpace , basis="canonical" )
    basis=="canonical" && (basis = CanonicalBasis( hilbert ))
    op = Operator( SymbolCreationOperator((snospin...,"u")) , basis )*
         Operator( SymbolAnnihilationOperator((snospin...,"d")) , basis )+
         Operator( SymbolCreationOperator((snospin...,"d")) , basis )*
         Operator( SymbolAnnihilationOperator((snospin...,"u")) , basis )
    return op 
end
spin_x( snospin::Tuple , hilbert::HilbertSpace , basis="canonical" ) = 
    0.5 * pauli_x( snospin , hilbert )

function pauli_y( snospin::Tuple , hilbert::HilbertSpace )
    basis=="canonical" && (basis = CanonicalBasis( hilbert ))
    op = Operator( SymbolCreationOperator((snospin...,"u")) , basis )*
         Operator( SymbolAnnihilationOperator((snospin...,"d")) , basis )-
         Operator( SymbolCreationOperator((snospin...,"d")) , basis )*
         Operator( SymbolAnnihilationOperator((snospin...,"u")) , basis )
    op *= -1.0im
    return op 
end
spin_y( snospin::Tuple , hilbert::HilbertSpace , basis="canonical" ) = 
    0.5 * pauli_y( snospin , hilbert , basis )

pauli2( snospin::Tuple , hilbert::HilbertSpace , basis="canonical" ) = 
                                                   pauli_x(snospin,hilbert)^2 + 
                                                   pauli_y(snospin,hilbert)^2 + 
                                                   pauli_z(snospin,hilbert)^2 
spin2( snospin::Tuple , hilbert::HilbertSpace , basis ) = 
    0.25 * pauli2( snospin , hilbert )

function pauli2_total( hilbert::HilbertSpace , basis="canonical" )
    snospins = Set([s[1:3] for s in hilbert.states])
    basis=="canonical" && (basis = CanonicalBasis( hilbert ))
    px = Operator( 0 , basis ) 
    py = Operator( 0 , basis )
    pz = Operator( 0 , basis )
    for snospin in snospins 
        px += pauli_x( snospin , hilbert , basis )
        py += pauli_y( snospin , hilbert , basis )
        pz += pauli_z( snospin , hilbert , basis )
    end
    return px^2 + py^2 + pz^2
end
spin2_total( hilbert::HilbertSpace , basis="canonical" ) =
    0.25 * pauli2_total( hilbert )


# ###### #
# NUMBER #
# ...... #
function number( stuple::Tuple , basis::CanonicalBasis )
    if length(stuple)==4
        c = Operator( SymbolCreationOperator(stuple) , basis )
        a = Operator( SymbolAnnihilationOperator(stuple) , basis )
        return Operator( SymbolCreationOperator(stuple) , basis )*
               Operator( SymbolAnnihilationOperator(stuple) , basis )
    elseif length(stuple)==3
        sup = ( stuple... , "u" )
        sdo = ( stuple... , "d" )
        #cup = Operator( SymbolCreationOperator(sup) , basis )
        #cdo = Operator( SymbolCreationOperator(sdo) , basis )
        #aup = Operator( SymbolAnnihilationOperator(sup) , basis )
        #ado = Operator( SymbolAnnihilationOperator(sdo) , basis )
        return Operator( SymbolCreationOperator(sup) , basis , herm=true )*
               Operator( SymbolAnnihilationOperator(sup) , basis , herm=true )+
               Operator( SymbolCreationOperator(sdo) , basis , herm=true )*
               Operator( SymbolAnnihilationOperator(sdo) , basis , herm=true )
        #return cup*aup + cdo*ado
    end
end
function number( hilbert::HilbertSpace )
    snospins = Set([s[1:3] for s in hilbert.states])
    basis = CanonicalBasis( hilbert )
    n = Operator( 0 , basis )
    for snospin in snospins
        n += number( snospin , basis )
    end
    return n 
end


# ####### #
# HOPPING #
# ....... #
function hop_1to2( snospin1 , snospin2 , basis )
    c_2u = Operator( SymbolCreationOperator((snospin2...,"u")) , basis )
    c_2d = Operator( SymbolCreationOperator((snospin2...,"d")) , basis )
    a_1u = Operator( SymbolAnnihilationOperator((snospin1...,"u")) , basis )
    a_1d = Operator( SymbolAnnihilationOperator((snospin1...,"d")) , basis )
    return c_2u*a_1u + c_2d*a_1d
end

function hop_bidirectional( snospin1 , snospin2 , basis )
    c_1u = Operator( SymbolCreationOperator((snospin1...,"u")) , basis )
    c_1d = Operator( SymbolCreationOperator((snospin1...,"d")) , basis )
    c_2u = Operator( SymbolCreationOperator((snospin2...,"u")) , basis )
    c_2d = Operator( SymbolCreationOperator((snospin2...,"d")) , basis )
    a_1u = Operator( SymbolAnnihilationOperator((snospin1...,"u")) , basis )
    a_1d = Operator( SymbolAnnihilationOperator((snospin1...,"d")) , basis )
    a_2u = Operator( SymbolAnnihilationOperator((snospin2...,"u")) , basis )
    a_2d = Operator( SymbolAnnihilationOperator((snospin2...,"d")) , basis )
    return c_1u*a_2u + c_1d*a_2d + c_2u*a_1u + c_2d*a_1d
end
