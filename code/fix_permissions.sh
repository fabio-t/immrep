find . -user $USER -type d -exec chmod 775 {} \; -exec chgrp bioinfo {} \; -exec chmod g+s {} \;
find . -user $USER -type f -exec chmod 664 {} \; -exec chgrp bioinfo {} \;

