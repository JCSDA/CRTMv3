MODULE CSEM_String_Utility
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: Str2UpCase
   PUBLIC :: Str2LowCase

   CHARACTER( * ), PRIVATE, PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
   CHARACTER( * ), PRIVATE, PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

CONTAINS

  FUNCTION Str2UpCase ( Input_String ) RESULT ( Output_String )
    ! -- Argument and result
    CHARACTER( * ), INTENT( IN )     :: Input_String
    CHARACTER( LEN( Input_String ) ) :: Output_String
    ! -- Local variables
    INTEGER :: i, n

    ! -- Copy input string
    Output_String = Input_String
    ! -- Loop over string elements
    DO i = 1, LEN( Output_String )
       ! -- Find location of letter in lower case constant string
       n = INDEX( LOWER_CASE, Output_String( i:i ) )
       ! -- If current substring is a lower case letter, make it upper case
       IF ( n /= 0 ) Output_String( i:i ) = UPPER_CASE( n:n )
    END DO
  END FUNCTION Str2UpCase

  FUNCTION Str2LowCase ( Input_String ) RESULT ( Output_String )
    ! -- Argument and result
    CHARACTER( * ), INTENT( IN )     :: Input_String
    CHARACTER( LEN( Input_String ) ) :: Output_String
    ! -- Local variables
    INTEGER :: i, n

    ! -- Copy input string
    Output_String = Input_String
    ! -- Loop over string elements
    DO i = 1, LEN( Output_String )
       ! -- Find location of letter in upper case constant string
       n = INDEX( UPPER_CASE, Output_String( i:i ) )
       ! -- If current substring is an upper case letter, make it lower case
       IF ( n /= 0 ) Output_String( i:i ) = LOWER_CASE( n:n )
    END DO
  END FUNCTION Str2LowCase
END MODULE CSEM_String_Utility 
